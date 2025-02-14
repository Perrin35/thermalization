OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.90647107) q[0];
sx q[0];
rz(4.8405092) q[0];
sx q[0];
rz(9.7066896) q[0];
rz(-2.6126722) q[1];
sx q[1];
rz(-1.6493874) q[1];
sx q[1];
rz(-1.5780916) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.603005) q[0];
sx q[0];
rz(-1.2055802) q[0];
sx q[0];
rz(0.19293789) q[0];
rz(-pi) q[1];
rz(2.8667502) q[2];
sx q[2];
rz(-1.5759528) q[2];
sx q[2];
rz(2.8066563) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.2927478) q[1];
sx q[1];
rz(-0.94460058) q[1];
sx q[1];
rz(1.7262154) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2658843) q[3];
sx q[3];
rz(-0.90581363) q[3];
sx q[3];
rz(0.18892442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6866744) q[2];
sx q[2];
rz(-2.8439971) q[2];
sx q[2];
rz(2.5285524) q[2];
rz(-0.4736627) q[3];
sx q[3];
rz(-1.9434171) q[3];
sx q[3];
rz(-1.4325498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4050201) q[0];
sx q[0];
rz(-0.18017811) q[0];
sx q[0];
rz(-0.75102425) q[0];
rz(-0.48149064) q[1];
sx q[1];
rz(-1.0844237) q[1];
sx q[1];
rz(-0.96985936) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6365627) q[0];
sx q[0];
rz(-1.9410994) q[0];
sx q[0];
rz(-2.4634393) q[0];
rz(2.9257751) q[2];
sx q[2];
rz(-1.6105798) q[2];
sx q[2];
rz(1.7173187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6578601) q[1];
sx q[1];
rz(-2.2071903) q[1];
sx q[1];
rz(2.8191393) q[1];
rz(-2.5327998) q[3];
sx q[3];
rz(-2.4498508) q[3];
sx q[3];
rz(2.5562654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91886175) q[2];
sx q[2];
rz(-1.6657882) q[2];
sx q[2];
rz(-2.6027423) q[2];
rz(-3.1239964) q[3];
sx q[3];
rz(-2.9774057) q[3];
sx q[3];
rz(-0.9084475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70055145) q[0];
sx q[0];
rz(-0.38527641) q[0];
sx q[0];
rz(0.061263099) q[0];
rz(1.6368658) q[1];
sx q[1];
rz(-2.442339) q[1];
sx q[1];
rz(1.1963074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.969097) q[0];
sx q[0];
rz(-0.75706702) q[0];
sx q[0];
rz(-1.2189381) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5012791) q[2];
sx q[2];
rz(-1.3090773) q[2];
sx q[2];
rz(1.6659322) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7972826) q[1];
sx q[1];
rz(-0.6629262) q[1];
sx q[1];
rz(1.0271038) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7964304) q[3];
sx q[3];
rz(-2.2335839) q[3];
sx q[3];
rz(1.8927285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4521744) q[2];
sx q[2];
rz(-1.3070561) q[2];
sx q[2];
rz(2.7090731) q[2];
rz(-1.0906667) q[3];
sx q[3];
rz(-0.64041036) q[3];
sx q[3];
rz(-2.045491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5525621) q[0];
sx q[0];
rz(-2.2046389) q[0];
sx q[0];
rz(0.20137782) q[0];
rz(-0.51678139) q[1];
sx q[1];
rz(-0.38025451) q[1];
sx q[1];
rz(-0.94863272) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1051335) q[0];
sx q[0];
rz(-2.3984342) q[0];
sx q[0];
rz(-2.0997597) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2791512) q[2];
sx q[2];
rz(-1.3459599) q[2];
sx q[2];
rz(0.98029691) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3310254) q[1];
sx q[1];
rz(-2.0750891) q[1];
sx q[1];
rz(-1.585998) q[1];
rz(-pi) q[2];
x q[2];
rz(1.473677) q[3];
sx q[3];
rz(-2.0896974) q[3];
sx q[3];
rz(-0.023954602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2739233) q[2];
sx q[2];
rz(-0.40931585) q[2];
sx q[2];
rz(-0.55740994) q[2];
rz(-3.0607767) q[3];
sx q[3];
rz(-1.2405453) q[3];
sx q[3];
rz(1.2062937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7252561) q[0];
sx q[0];
rz(-1.6660322) q[0];
sx q[0];
rz(-2.1112554) q[0];
rz(-2.985785) q[1];
sx q[1];
rz(-1.067433) q[1];
sx q[1];
rz(-2.9373998) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1847398) q[0];
sx q[0];
rz(-2.889688) q[0];
sx q[0];
rz(2.8624318) q[0];
rz(0.5459909) q[2];
sx q[2];
rz(-0.24790774) q[2];
sx q[2];
rz(1.966983) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5203917) q[1];
sx q[1];
rz(-1.0869496) q[1];
sx q[1];
rz(1.8041759) q[1];
x q[2];
rz(0.24498265) q[3];
sx q[3];
rz(-0.7881654) q[3];
sx q[3];
rz(-0.80313166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2121409) q[2];
sx q[2];
rz(-1.6484304) q[2];
sx q[2];
rz(0.41352752) q[2];
rz(-2.2996969) q[3];
sx q[3];
rz(-1.7588408) q[3];
sx q[3];
rz(-0.97257417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0464762) q[0];
sx q[0];
rz(-1.0914509) q[0];
sx q[0];
rz(0.0055775642) q[0];
rz(1.3439517) q[1];
sx q[1];
rz(-2.9426212) q[1];
sx q[1];
rz(2.4406348) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.287375) q[0];
sx q[0];
rz(-1.1297884) q[0];
sx q[0];
rz(0.99325755) q[0];
x q[1];
rz(-2.2118819) q[2];
sx q[2];
rz(-1.9308917) q[2];
sx q[2];
rz(1.4907229) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4206767) q[1];
sx q[1];
rz(-1.5297959) q[1];
sx q[1];
rz(0.27537217) q[1];
rz(-pi) q[2];
rz(-2.5684729) q[3];
sx q[3];
rz(-1.1240715) q[3];
sx q[3];
rz(-0.61441159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.31412101) q[2];
sx q[2];
rz(-0.33385971) q[2];
sx q[2];
rz(1.9742879) q[2];
rz(-2.2589034) q[3];
sx q[3];
rz(-0.76789951) q[3];
sx q[3];
rz(-0.32514969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2348787) q[0];
sx q[0];
rz(-1.119708) q[0];
sx q[0];
rz(0.085302189) q[0];
rz(-0.75434297) q[1];
sx q[1];
rz(-2.6383548) q[1];
sx q[1];
rz(-2.1139961) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2404382) q[0];
sx q[0];
rz(-1.580997) q[0];
sx q[0];
rz(2.9290694) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1336537) q[2];
sx q[2];
rz(-2.5449341) q[2];
sx q[2];
rz(-0.82574979) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4484549) q[1];
sx q[1];
rz(-1.692728) q[1];
sx q[1];
rz(1.8961402) q[1];
rz(-pi) q[2];
rz(1.0124765) q[3];
sx q[3];
rz(-1.9838247) q[3];
sx q[3];
rz(0.2193887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9895642) q[2];
sx q[2];
rz(-0.98735183) q[2];
sx q[2];
rz(2.7222471) q[2];
rz(-0.63465214) q[3];
sx q[3];
rz(-0.80497634) q[3];
sx q[3];
rz(1.1254719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80970508) q[0];
sx q[0];
rz(-2.2655847) q[0];
sx q[0];
rz(-0.69212717) q[0];
rz(0.57811111) q[1];
sx q[1];
rz(-2.7222241) q[1];
sx q[1];
rz(-1.3828166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0771187) q[0];
sx q[0];
rz(-0.070264272) q[0];
sx q[0];
rz(1.8514567) q[0];
rz(-pi) q[1];
rz(0.67202576) q[2];
sx q[2];
rz(-2.1161072) q[2];
sx q[2];
rz(-1.3915075) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9005712) q[1];
sx q[1];
rz(-1.3374778) q[1];
sx q[1];
rz(2.1395348) q[1];
x q[2];
rz(2.0973849) q[3];
sx q[3];
rz(-2.3447737) q[3];
sx q[3];
rz(-0.45180333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.77070037) q[2];
sx q[2];
rz(-0.53052491) q[2];
sx q[2];
rz(1.4250866) q[2];
rz(-0.51618451) q[3];
sx q[3];
rz(-2.1631212) q[3];
sx q[3];
rz(-1.9019351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4952963) q[0];
sx q[0];
rz(-1.7072059) q[0];
sx q[0];
rz(0.37149757) q[0];
rz(2.3067572) q[1];
sx q[1];
rz(-2.3284262) q[1];
sx q[1];
rz(3.1239948) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8397578) q[0];
sx q[0];
rz(-2.5389414) q[0];
sx q[0];
rz(1.0687123) q[0];
x q[1];
rz(0.27522343) q[2];
sx q[2];
rz(-0.75481269) q[2];
sx q[2];
rz(1.7856904) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0068215) q[1];
sx q[1];
rz(-2.4398514) q[1];
sx q[1];
rz(-0.59533466) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60815717) q[3];
sx q[3];
rz(-1.8377234) q[3];
sx q[3];
rz(1.9723231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1060433) q[2];
sx q[2];
rz(-0.86842662) q[2];
sx q[2];
rz(3.0412728) q[2];
rz(-0.22942461) q[3];
sx q[3];
rz(-1.0474297) q[3];
sx q[3];
rz(-0.44982287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4556274) q[0];
sx q[0];
rz(-0.28814155) q[0];
sx q[0];
rz(-2.567754) q[0];
rz(2.0876743) q[1];
sx q[1];
rz(-1.046448) q[1];
sx q[1];
rz(0.67646772) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5887816) q[0];
sx q[0];
rz(-1.2924606) q[0];
sx q[0];
rz(1.3180483) q[0];
x q[1];
rz(-1.1660812) q[2];
sx q[2];
rz(-1.8076483) q[2];
sx q[2];
rz(2.7450457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0770196) q[1];
sx q[1];
rz(-1.4823682) q[1];
sx q[1];
rz(2.6131373) q[1];
rz(-2.1031688) q[3];
sx q[3];
rz(-1.3992953) q[3];
sx q[3];
rz(0.72491437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.68710589) q[2];
sx q[2];
rz(-0.7578308) q[2];
sx q[2];
rz(-1.494361) q[2];
rz(-2.2729661) q[3];
sx q[3];
rz(-0.70610154) q[3];
sx q[3];
rz(-0.47006616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0710707) q[0];
sx q[0];
rz(-1.1916397) q[0];
sx q[0];
rz(1.2217039) q[0];
rz(1.8078049) q[1];
sx q[1];
rz(-1.8232657) q[1];
sx q[1];
rz(2.5008536) q[1];
rz(2.376198) q[2];
sx q[2];
rz(-1.6366048) q[2];
sx q[2];
rz(-0.66726782) q[2];
rz(-1.3258237) q[3];
sx q[3];
rz(-0.180937) q[3];
sx q[3];
rz(-2.141249) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
