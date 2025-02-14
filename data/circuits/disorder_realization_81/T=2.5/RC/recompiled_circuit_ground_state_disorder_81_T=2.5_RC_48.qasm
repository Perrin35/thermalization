OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9658907) q[0];
sx q[0];
rz(-1.8718636) q[0];
sx q[0];
rz(1.1748535) q[0];
rz(-1.084561) q[1];
sx q[1];
rz(-1.4718055) q[1];
sx q[1];
rz(0.56301277) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36328735) q[0];
sx q[0];
rz(-1.7666923) q[0];
sx q[0];
rz(-1.6250074) q[0];
x q[1];
rz(2.1839574) q[2];
sx q[2];
rz(-2.4102825) q[2];
sx q[2];
rz(2.4201408) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1606428) q[1];
sx q[1];
rz(-1.8997571) q[1];
sx q[1];
rz(-1.739722) q[1];
x q[2];
rz(0.38083846) q[3];
sx q[3];
rz(-1.2361174) q[3];
sx q[3];
rz(-0.0045485529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53434831) q[2];
sx q[2];
rz(-1.6515942) q[2];
sx q[2];
rz(1.6687757) q[2];
rz(1.8913174) q[3];
sx q[3];
rz(-1.5273124) q[3];
sx q[3];
rz(0.97243029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2431353) q[0];
sx q[0];
rz(-0.055483015) q[0];
sx q[0];
rz(-2.6836416) q[0];
rz(3.0398439) q[1];
sx q[1];
rz(-1.190217) q[1];
sx q[1];
rz(2.8153458) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3889623) q[0];
sx q[0];
rz(-1.2241619) q[0];
sx q[0];
rz(2.017157) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8433766) q[2];
sx q[2];
rz(-1.3591849) q[2];
sx q[2];
rz(-0.22842184) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5811903) q[1];
sx q[1];
rz(-0.96424864) q[1];
sx q[1];
rz(0.62967695) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74298782) q[3];
sx q[3];
rz(-1.3338193) q[3];
sx q[3];
rz(-1.5239511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0094256224) q[2];
sx q[2];
rz(-2.1475809) q[2];
sx q[2];
rz(1.2635292) q[2];
rz(3.1215014) q[3];
sx q[3];
rz(-0.76954904) q[3];
sx q[3];
rz(1.2300389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024071368) q[0];
sx q[0];
rz(-3.0776403) q[0];
sx q[0];
rz(-0.57620311) q[0];
rz(-1.6974712) q[1];
sx q[1];
rz(-1.2497808) q[1];
sx q[1];
rz(1.8796657) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5844628) q[0];
sx q[0];
rz(-0.31079159) q[0];
sx q[0];
rz(-1.6755098) q[0];
rz(-0.70749486) q[2];
sx q[2];
rz(-0.65366483) q[2];
sx q[2];
rz(-2.737239) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60186767) q[1];
sx q[1];
rz(-2.4721708) q[1];
sx q[1];
rz(0.4008534) q[1];
rz(1.0141054) q[3];
sx q[3];
rz(-2.2051349) q[3];
sx q[3];
rz(1.9012828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5835517) q[2];
sx q[2];
rz(-0.59528196) q[2];
sx q[2];
rz(-0.15312791) q[2];
rz(0.88614744) q[3];
sx q[3];
rz(-1.359442) q[3];
sx q[3];
rz(-1.9396293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8837638) q[0];
sx q[0];
rz(-1.184329) q[0];
sx q[0];
rz(-2.9828239) q[0];
rz(2.1067045) q[1];
sx q[1];
rz(-2.1885469) q[1];
sx q[1];
rz(0.75611702) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6203228) q[0];
sx q[0];
rz(-1.5128994) q[0];
sx q[0];
rz(-1.6211364) q[0];
rz(1.6768084) q[2];
sx q[2];
rz(-0.25463018) q[2];
sx q[2];
rz(0.98843304) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9489574) q[1];
sx q[1];
rz(-1.3657943) q[1];
sx q[1];
rz(2.8947049) q[1];
rz(-pi) q[2];
rz(2.4231195) q[3];
sx q[3];
rz(-0.30820307) q[3];
sx q[3];
rz(-2.001345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5990344) q[2];
sx q[2];
rz(-0.90347806) q[2];
sx q[2];
rz(0.8503882) q[2];
rz(-1.4706069) q[3];
sx q[3];
rz(-1.8149523) q[3];
sx q[3];
rz(-2.3175122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65664148) q[0];
sx q[0];
rz(-1.739772) q[0];
sx q[0];
rz(2.9630419) q[0];
rz(1.8932331) q[1];
sx q[1];
rz(-1.6105885) q[1];
sx q[1];
rz(0.53370968) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8781136) q[0];
sx q[0];
rz(-1.7537259) q[0];
sx q[0];
rz(0.35533743) q[0];
rz(-pi) q[1];
rz(2.2703993) q[2];
sx q[2];
rz(-1.9661926) q[2];
sx q[2];
rz(-1.7311387) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43161422) q[1];
sx q[1];
rz(-1.5182765) q[1];
sx q[1];
rz(-2.2389328) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7469095) q[3];
sx q[3];
rz(-0.51100547) q[3];
sx q[3];
rz(0.7738302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8278213) q[2];
sx q[2];
rz(-1.2510108) q[2];
sx q[2];
rz(-2.0884464) q[2];
rz(0.18038067) q[3];
sx q[3];
rz(-0.80612055) q[3];
sx q[3];
rz(-0.36261305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5617705) q[0];
sx q[0];
rz(-1.0012015) q[0];
sx q[0];
rz(-1.2286105) q[0];
rz(-0.46663943) q[1];
sx q[1];
rz(-1.2184315) q[1];
sx q[1];
rz(3.012588) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6138959) q[0];
sx q[0];
rz(-1.6718699) q[0];
sx q[0];
rz(1.4368426) q[0];
x q[1];
rz(-2.9797709) q[2];
sx q[2];
rz(-2.2941049) q[2];
sx q[2];
rz(2.3486193) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58632233) q[1];
sx q[1];
rz(-0.96809371) q[1];
sx q[1];
rz(-1.210825) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3025593) q[3];
sx q[3];
rz(-1.7058289) q[3];
sx q[3];
rz(2.6756833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1351607) q[2];
sx q[2];
rz(-2.5401523) q[2];
sx q[2];
rz(-0.50430164) q[2];
rz(-2.048061) q[3];
sx q[3];
rz(-1.8176327) q[3];
sx q[3];
rz(0.95660153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4920376) q[0];
sx q[0];
rz(-1.858359) q[0];
sx q[0];
rz(1.5851703) q[0];
rz(-0.57836142) q[1];
sx q[1];
rz(-1.2588986) q[1];
sx q[1];
rz(-0.10231054) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1083173) q[0];
sx q[0];
rz(-0.96848291) q[0];
sx q[0];
rz(2.6541416) q[0];
rz(-pi) q[1];
rz(-0.39686508) q[2];
sx q[2];
rz(-1.3584653) q[2];
sx q[2];
rz(-3.120829) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31033406) q[1];
sx q[1];
rz(-1.4142904) q[1];
sx q[1];
rz(-2.8328647) q[1];
rz(-pi) q[2];
rz(2.3897947) q[3];
sx q[3];
rz(-2.2357142) q[3];
sx q[3];
rz(-0.81937688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9461225) q[2];
sx q[2];
rz(-1.8778233) q[2];
sx q[2];
rz(-2.5275687) q[2];
rz(-0.92769235) q[3];
sx q[3];
rz(-1.5427019) q[3];
sx q[3];
rz(-0.6774925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55352655) q[0];
sx q[0];
rz(-0.29356846) q[0];
sx q[0];
rz(-1.9422096) q[0];
rz(-1.9623494) q[1];
sx q[1];
rz(-2.113138) q[1];
sx q[1];
rz(0.95338043) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48483585) q[0];
sx q[0];
rz(-1.0547045) q[0];
sx q[0];
rz(2.1950095) q[0];
rz(-1.751975) q[2];
sx q[2];
rz(-2.5963514) q[2];
sx q[2];
rz(-0.94961005) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7443707) q[1];
sx q[1];
rz(-0.83852856) q[1];
sx q[1];
rz(-2.9930755) q[1];
rz(1.6497506) q[3];
sx q[3];
rz(-1.0903666) q[3];
sx q[3];
rz(1.3955658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.392268) q[2];
sx q[2];
rz(-3.0159123) q[2];
sx q[2];
rz(-0.18088642) q[2];
rz(1.1130029) q[3];
sx q[3];
rz(-1.0180611) q[3];
sx q[3];
rz(-0.38016144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0893294) q[0];
sx q[0];
rz(-0.90173975) q[0];
sx q[0];
rz(-2.3774636) q[0];
rz(-2.1185421) q[1];
sx q[1];
rz(-1.6108395) q[1];
sx q[1];
rz(1.4952362) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8500124) q[0];
sx q[0];
rz(-0.33667013) q[0];
sx q[0];
rz(2.6597937) q[0];
x q[1];
rz(-2.1713397) q[2];
sx q[2];
rz(-1.8489328) q[2];
sx q[2];
rz(-2.8153265) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0937522) q[1];
sx q[1];
rz(-1.6421659) q[1];
sx q[1];
rz(2.755295) q[1];
x q[2];
rz(0.44166724) q[3];
sx q[3];
rz(-1.212647) q[3];
sx q[3];
rz(1.4513858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9493147) q[2];
sx q[2];
rz(-1.01769) q[2];
sx q[2];
rz(-2.6940572) q[2];
rz(3.0053075) q[3];
sx q[3];
rz(-0.49301967) q[3];
sx q[3];
rz(-0.59806699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3560781) q[0];
sx q[0];
rz(-1.0250174) q[0];
sx q[0];
rz(-2.6688975) q[0];
rz(0.39974943) q[1];
sx q[1];
rz(-0.70416299) q[1];
sx q[1];
rz(-2.3230816) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2855354) q[0];
sx q[0];
rz(-1.7489079) q[0];
sx q[0];
rz(2.8277603) q[0];
rz(1.1498951) q[2];
sx q[2];
rz(-0.93846655) q[2];
sx q[2];
rz(2.9078751) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1295197) q[1];
sx q[1];
rz(-2.8188779) q[1];
sx q[1];
rz(1.580915) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7943939) q[3];
sx q[3];
rz(-1.7337017) q[3];
sx q[3];
rz(2.3293975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5588536) q[2];
sx q[2];
rz(-1.0820729) q[2];
sx q[2];
rz(1.2782798) q[2];
rz(1.9418955) q[3];
sx q[3];
rz(-1.4408828) q[3];
sx q[3];
rz(0.27967683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3661135) q[0];
sx q[0];
rz(-1.4237325) q[0];
sx q[0];
rz(1.7970418) q[0];
rz(-1.7300425) q[1];
sx q[1];
rz(-2.330214) q[1];
sx q[1];
rz(-0.65705962) q[1];
rz(-3.0006164) q[2];
sx q[2];
rz(-1.5700414) q[2];
sx q[2];
rz(1.0667917) q[2];
rz(2.2349002) q[3];
sx q[3];
rz(-1.8421052) q[3];
sx q[3];
rz(2.6143034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
