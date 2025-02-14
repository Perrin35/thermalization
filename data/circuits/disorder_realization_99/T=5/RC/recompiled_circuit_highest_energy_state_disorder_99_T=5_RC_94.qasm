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
rz(-1.7412269) q[0];
sx q[0];
rz(-2.9018612) q[0];
sx q[0];
rz(2.8170407) q[0];
rz(2.1895154) q[1];
sx q[1];
rz(-2.9157186) q[1];
sx q[1];
rz(-1.3083375) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4885534) q[0];
sx q[0];
rz(-1.4346315) q[0];
sx q[0];
rz(-0.67739886) q[0];
rz(-pi) q[1];
rz(0.24271528) q[2];
sx q[2];
rz(-2.6449892) q[2];
sx q[2];
rz(2.9951823) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5950299) q[1];
sx q[1];
rz(-2.0339222) q[1];
sx q[1];
rz(-0.28966499) q[1];
rz(-2.4423787) q[3];
sx q[3];
rz(-1.6788071) q[3];
sx q[3];
rz(2.8349769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41210458) q[2];
sx q[2];
rz(-1.8821913) q[2];
sx q[2];
rz(-2.7903902) q[2];
rz(1.2900194) q[3];
sx q[3];
rz(-2.0100287) q[3];
sx q[3];
rz(-1.9180408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7555162) q[0];
sx q[0];
rz(-1.2048683) q[0];
sx q[0];
rz(0.93038857) q[0];
rz(-1.3049841) q[1];
sx q[1];
rz(-0.84140673) q[1];
sx q[1];
rz(2.5321541) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7284645) q[0];
sx q[0];
rz(-0.43914686) q[0];
sx q[0];
rz(1.3261611) q[0];
rz(-pi) q[1];
rz(-1.2728884) q[2];
sx q[2];
rz(-1.6570083) q[2];
sx q[2];
rz(0.49315587) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8011368) q[1];
sx q[1];
rz(-1.1994455) q[1];
sx q[1];
rz(0.39210971) q[1];
x q[2];
rz(1.3136765) q[3];
sx q[3];
rz(-2.5559396) q[3];
sx q[3];
rz(0.7440834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.096752) q[2];
sx q[2];
rz(-2.1614306) q[2];
sx q[2];
rz(0.074782221) q[2];
rz(-0.52037248) q[3];
sx q[3];
rz(-2.4986391) q[3];
sx q[3];
rz(-0.61409942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5388913) q[0];
sx q[0];
rz(-2.4245872) q[0];
sx q[0];
rz(0.30211788) q[0];
rz(2.865454) q[1];
sx q[1];
rz(-1.3172251) q[1];
sx q[1];
rz(-1.709323) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30423588) q[0];
sx q[0];
rz(-0.68528995) q[0];
sx q[0];
rz(0.94288148) q[0];
x q[1];
rz(2.7652241) q[2];
sx q[2];
rz(-0.62877765) q[2];
sx q[2];
rz(-0.16964682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.773743) q[1];
sx q[1];
rz(-2.5605695) q[1];
sx q[1];
rz(-1.2817205) q[1];
x q[2];
rz(2.0366482) q[3];
sx q[3];
rz(-1.5072952) q[3];
sx q[3];
rz(2.7257763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8351195) q[2];
sx q[2];
rz(-0.46572954) q[2];
sx q[2];
rz(1.845537) q[2];
rz(0.45201388) q[3];
sx q[3];
rz(-1.1072423) q[3];
sx q[3];
rz(0.38255102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52962676) q[0];
sx q[0];
rz(-1.9744248) q[0];
sx q[0];
rz(-1.3324598) q[0];
rz(-2.0924856) q[1];
sx q[1];
rz(-2.0499947) q[1];
sx q[1];
rz(1.5886935) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6352672) q[0];
sx q[0];
rz(-2.679279) q[0];
sx q[0];
rz(-2.4145831) q[0];
rz(-pi) q[1];
x q[1];
rz(1.595644) q[2];
sx q[2];
rz(-2.6406807) q[2];
sx q[2];
rz(2.3786169) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3192352) q[1];
sx q[1];
rz(-2.2062782) q[1];
sx q[1];
rz(0.24680071) q[1];
rz(-pi) q[2];
rz(2.6060206) q[3];
sx q[3];
rz(-0.93129292) q[3];
sx q[3];
rz(0.65062338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66854746) q[2];
sx q[2];
rz(-2.1210402) q[2];
sx q[2];
rz(-0.2505396) q[2];
rz(-0.9969095) q[3];
sx q[3];
rz(-1.1397811) q[3];
sx q[3];
rz(-0.38058773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8119891) q[0];
sx q[0];
rz(-2.6530837) q[0];
sx q[0];
rz(1.3478152) q[0];
rz(0.40428058) q[1];
sx q[1];
rz(-0.75899044) q[1];
sx q[1];
rz(-1.1291198) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4816138) q[0];
sx q[0];
rz(-2.3140175) q[0];
sx q[0];
rz(-1.0218234) q[0];
x q[1];
rz(1.2329699) q[2];
sx q[2];
rz(-2.4475636) q[2];
sx q[2];
rz(2.239486) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.68539219) q[1];
sx q[1];
rz(-0.90185748) q[1];
sx q[1];
rz(-3.0078956) q[1];
rz(-0.31480883) q[3];
sx q[3];
rz(-1.4658576) q[3];
sx q[3];
rz(0.85280692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79444844) q[2];
sx q[2];
rz(-0.930154) q[2];
sx q[2];
rz(1.8974737) q[2];
rz(2.0813023) q[3];
sx q[3];
rz(-2.5131707) q[3];
sx q[3];
rz(-0.30317831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65619549) q[0];
sx q[0];
rz(-3.0363016) q[0];
sx q[0];
rz(2.7144077) q[0];
rz(3.1071013) q[1];
sx q[1];
rz(-1.5723615) q[1];
sx q[1];
rz(3.131386) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0611864) q[0];
sx q[0];
rz(-1.5725333) q[0];
sx q[0];
rz(-3.1392359) q[0];
rz(-pi) q[1];
rz(-2.9811007) q[2];
sx q[2];
rz(-2.1969323) q[2];
sx q[2];
rz(0.4753091) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3349377) q[1];
sx q[1];
rz(-0.68938556) q[1];
sx q[1];
rz(1.4695808) q[1];
rz(-2.4236492) q[3];
sx q[3];
rz(-2.4417851) q[3];
sx q[3];
rz(-1.0216624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61937845) q[2];
sx q[2];
rz(-0.41316119) q[2];
sx q[2];
rz(-2.348032) q[2];
rz(0.86269745) q[3];
sx q[3];
rz(-1.5669275) q[3];
sx q[3];
rz(2.4376552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70235395) q[0];
sx q[0];
rz(-0.84504253) q[0];
sx q[0];
rz(-2.0080361) q[0];
rz(-2.0639065) q[1];
sx q[1];
rz(-2.6262941) q[1];
sx q[1];
rz(-0.30212197) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0724807) q[0];
sx q[0];
rz(-1.9078507) q[0];
sx q[0];
rz(1.891579) q[0];
x q[1];
rz(2.5100915) q[2];
sx q[2];
rz(-0.69984791) q[2];
sx q[2];
rz(1.8818784) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4430284) q[1];
sx q[1];
rz(-2.5686567) q[1];
sx q[1];
rz(-0.40614508) q[1];
rz(-pi) q[2];
rz(-2.1572635) q[3];
sx q[3];
rz(-0.97198379) q[3];
sx q[3];
rz(2.5564155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2332396) q[2];
sx q[2];
rz(-2.2519604) q[2];
sx q[2];
rz(1.3471777) q[2];
rz(-1.8917278) q[3];
sx q[3];
rz(-1.9439387) q[3];
sx q[3];
rz(-3.0344322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5477448) q[0];
sx q[0];
rz(-2.7084454) q[0];
sx q[0];
rz(-0.10840848) q[0];
rz(1.0477061) q[1];
sx q[1];
rz(-2.3240418) q[1];
sx q[1];
rz(-1.0188867) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7664095) q[0];
sx q[0];
rz(-2.0926704) q[0];
sx q[0];
rz(0.96065582) q[0];
rz(2.3802929) q[2];
sx q[2];
rz(-1.5754964) q[2];
sx q[2];
rz(2.0705303) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0483886) q[1];
sx q[1];
rz(-2.0282413) q[1];
sx q[1];
rz(-2.0162575) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74014981) q[3];
sx q[3];
rz(-1.6751846) q[3];
sx q[3];
rz(-1.1889282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.63742796) q[2];
sx q[2];
rz(-0.97680682) q[2];
sx q[2];
rz(-0.43761474) q[2];
rz(0.49312433) q[3];
sx q[3];
rz(-2.2212641) q[3];
sx q[3];
rz(1.5639308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86579943) q[0];
sx q[0];
rz(-1.9075305) q[0];
sx q[0];
rz(0.27780521) q[0];
rz(1.5996784) q[1];
sx q[1];
rz(-2.1733687) q[1];
sx q[1];
rz(-0.42617282) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23754643) q[0];
sx q[0];
rz(-1.5874366) q[0];
sx q[0];
rz(-0.017701935) q[0];
rz(-pi) q[1];
rz(-2.2913646) q[2];
sx q[2];
rz(-0.93932187) q[2];
sx q[2];
rz(-0.15635083) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8746169) q[1];
sx q[1];
rz(-0.24493453) q[1];
sx q[1];
rz(1.0670948) q[1];
rz(-pi) q[2];
rz(-2.3087193) q[3];
sx q[3];
rz(-0.98372059) q[3];
sx q[3];
rz(-0.47490109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6466732) q[2];
sx q[2];
rz(-0.16725954) q[2];
sx q[2];
rz(-2.0885928) q[2];
rz(0.35887512) q[3];
sx q[3];
rz(-1.5263298) q[3];
sx q[3];
rz(-1.0927965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(2.7786355) q[0];
sx q[0];
rz(-0.21720049) q[0];
sx q[0];
rz(-2.2208075) q[0];
rz(2.6329363) q[1];
sx q[1];
rz(-2.3743036) q[1];
sx q[1];
rz(-2.7210534) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79179278) q[0];
sx q[0];
rz(-1.0726439) q[0];
sx q[0];
rz(0.80043261) q[0];
x q[1];
rz(2.8996572) q[2];
sx q[2];
rz(-2.1872949) q[2];
sx q[2];
rz(0.37504196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9892726) q[1];
sx q[1];
rz(-0.97022431) q[1];
sx q[1];
rz(2.5737073) q[1];
rz(0.8798265) q[3];
sx q[3];
rz(-1.4925957) q[3];
sx q[3];
rz(-1.2134589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4113808) q[2];
sx q[2];
rz(-2.3513942) q[2];
sx q[2];
rz(-0.31663695) q[2];
rz(2.5824879) q[3];
sx q[3];
rz(-2.5059301) q[3];
sx q[3];
rz(-2.3234698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.485514) q[0];
sx q[0];
rz(-1.0777127) q[0];
sx q[0];
rz(1.0832473) q[0];
rz(2.6651233) q[1];
sx q[1];
rz(-2.10119) q[1];
sx q[1];
rz(1.4269921) q[1];
rz(-0.086924788) q[2];
sx q[2];
rz(-0.82155052) q[2];
sx q[2];
rz(-0.62053298) q[2];
rz(1.1875759) q[3];
sx q[3];
rz(-1.4837617) q[3];
sx q[3];
rz(-1.0401534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
