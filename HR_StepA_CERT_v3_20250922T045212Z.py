# HR Step A — v3: integrated zero-table reader + VK constants loader/suggester
# Usage:
#   - Set CONFIG["zeros_file"] to your certified zero list up to T0 (format: "idx gamma" or "gamma").
#   - Set CONFIG["C_R"] to the certified kernel remainder constant.
#   - Either set CONFIG["vk_constants_json"] with {"B_VK":..., "b_VK":..., "x1":...}
#     or let the script suggest NON-CERTIFIED values from the VK region constant (55.241 by default).
#
# Output: a JSON and a PDF with (C_bajo, C_alto, C_empalme, C_tot) and diagnostics.
#
# NOTE: The suggestion for (B_VK, b_VK) is *not* a proof; for certification, provide vk_constants_json.
