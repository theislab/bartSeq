import copy
import math
import random
from logging import getLogger
from typing import List, Tuple


log = getLogger(__name__)

Arrangement = Tuple[float, List[int], List[int]]


class Optimizer:
    def __init__(self, config, sequence_set):
        self.config = config
        self.sequence_set = sequence_set

    def f(self, selected_pairs, selected_amplicons) -> float:
        sum_mfes = 0
        for seq1, pair1 in enumerate(selected_pairs):
            for seq2, pair2 in enumerate(selected_pairs):
                sum_mfes += self.sequence_set[seq1].mfes[selected_amplicons[seq1]][
                    pair1
                ][seq2][selected_amplicons[seq2]][pair2]
        return sum_mfes / 2

    def optimize(self) -> List[Arrangement]:
        log.info("optimizing")
        max_ind = self.config.opt_steps
        max_temperature = min(self.config.opt_max_temp, self.config.opt_steps)

        amplicon_lengths: List[int] = []
        set_lengths: List[List[int]] = []

        for i, gene in enumerate(self.sequence_set):
            aset = gene.amplicons
            amplicon_lengths.append(len(aset))
            set_lengths.append([])
            for amplicon in aset:
                set_lengths[i].append(len(amplicon.primer_set))

        combinations: List[Arrangement] = []
        act_temperature = max_temperature

        w = []
        for i, gene in enumerate(self.sequence_set):
            w.append(random.randint(0, amplicon_lengths[i] - 1))

        v = []
        for i, seq in enumerate(self.sequence_set):
            v.append(random.randint(0, set_lengths[i][w[i]] - 1))

        temp_steps = math.floor(max_ind / (max_temperature + 1))

        no_change = 0

        i = 0
        while no_change < 1500:
            # if i % math.floor(max_ind/10) == 0:
            #     print ((i/math.floor(max_ind))*100,"%")

            if i % 1000 == 0:
                log.info(
                    "No change: %s, temp: %s, score(v, w) = %s",
                    no_change,
                    act_temperature,
                    self.f(v, w),
                )

            if i % temp_steps == 0 and act_temperature != 0:
                act_temperature += -1

            j = random.randint(0, len(v) - 1)
            changed = False

            # in 20% of the cases change the amplicon
            if random.random() < 0.2:
                for k in random.sample(
                    range(0, amplicon_lengths[j]), amplicon_lengths[j]
                ):
                    v_temp = copy.copy(v)
                    for l in random.sample(
                        range(0, set_lengths[j][k]), set_lengths[j][k]
                    ):
                        v_temp[j] = l
                        w_temp = copy.copy(w)
                        w_temp[j] = k

                        old_score = self.f(v, w)
                        new_score = self.f(v_temp, w_temp)

                        if new_score <= old_score:
                            w[j] = k
                            v[j] = l
                        elif act_temperature > 0 and math.exp(
                            (old_score - new_score) / act_temperature
                        ) > random.uniform(0, 1):
                            w[j] = k
                            v[j] = l

                        if new_score < old_score:
                            changed = True
            else:
                for k in random.sample(
                    range(0, set_lengths[j][w[j]]), set_lengths[j][w[j]]
                ):
                    v_temp = copy.copy(v)
                    v_temp[j] = k

                    old_score = self.f(v, w)
                    new_score = self.f(v_temp, w)

                    if new_score <= old_score:
                        v[j] = k
                    elif act_temperature > 0 and math.exp(
                        (old_score - new_score) / act_temperature
                    ) > random.uniform(0, 1):
                        v[j] = k

                    if new_score < old_score:
                        changed = True

            if changed:
                no_change = 0
            else:
                no_change += 1
            combinations.append((self.f(v, w), copy.copy(v), copy.copy(w)))
            i += 1

        unique_arrangements: List[Arrangement] = []
        arrangements = sorted(combinations, key=lambda x: x[0])

        for run in arrangements:
            if run not in unique_arrangements:
                unique_arrangements.append(run)

        return unique_arrangements
