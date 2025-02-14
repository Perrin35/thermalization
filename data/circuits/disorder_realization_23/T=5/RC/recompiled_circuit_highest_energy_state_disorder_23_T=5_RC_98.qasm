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
rz(0.55573207) q[0];
sx q[0];
rz(-1.8620123) q[0];
sx q[0];
rz(-0.32787856) q[0];
rz(-2.9887587) q[1];
sx q[1];
rz(-2.6522377) q[1];
sx q[1];
rz(2.1305003) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0965138) q[0];
sx q[0];
rz(-0.9849087) q[0];
sx q[0];
rz(-1.3784598) q[0];
rz(-pi) q[1];
rz(1.7494124) q[2];
sx q[2];
rz(-1.9733323) q[2];
sx q[2];
rz(2.2312763) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7465621) q[1];
sx q[1];
rz(-0.85828188) q[1];
sx q[1];
rz(1.0814352) q[1];
x q[2];
rz(-2.070363) q[3];
sx q[3];
rz(-0.46854436) q[3];
sx q[3];
rz(0.51817453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8622387) q[2];
sx q[2];
rz(-0.78247672) q[2];
sx q[2];
rz(1.5581101) q[2];
rz(-0.33485788) q[3];
sx q[3];
rz(-2.0575276) q[3];
sx q[3];
rz(0.28086942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8849477) q[0];
sx q[0];
rz(-2.5769951) q[0];
sx q[0];
rz(-2.8096492) q[0];
rz(2.7815107) q[1];
sx q[1];
rz(-1.8372476) q[1];
sx q[1];
rz(0.28775451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.677414) q[0];
sx q[0];
rz(-0.3218284) q[0];
sx q[0];
rz(-2.4610956) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7878868) q[2];
sx q[2];
rz(-0.76375095) q[2];
sx q[2];
rz(-1.1909831) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9336492) q[1];
sx q[1];
rz(-0.48743409) q[1];
sx q[1];
rz(-0.80923621) q[1];
rz(-pi) q[2];
rz(-0.61057456) q[3];
sx q[3];
rz(-2.861851) q[3];
sx q[3];
rz(2.4795462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.45431367) q[2];
sx q[2];
rz(-1.6087029) q[2];
sx q[2];
rz(-1.7506556) q[2];
rz(2.2823997) q[3];
sx q[3];
rz(-1.8682559) q[3];
sx q[3];
rz(2.9505762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.745382) q[0];
sx q[0];
rz(-1.6062382) q[0];
sx q[0];
rz(-0.23042738) q[0];
rz(-0.51741171) q[1];
sx q[1];
rz(-2.0473862) q[1];
sx q[1];
rz(0.28883019) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5342425) q[0];
sx q[0];
rz(-1.6301883) q[0];
sx q[0];
rz(-1.1618105) q[0];
rz(-0.46858139) q[2];
sx q[2];
rz(-1.6372674) q[2];
sx q[2];
rz(2.2867212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3687262) q[1];
sx q[1];
rz(-1.9802367) q[1];
sx q[1];
rz(-0.057842908) q[1];
x q[2];
rz(-2.9656319) q[3];
sx q[3];
rz(-1.3506512) q[3];
sx q[3];
rz(-2.1484321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0383703) q[2];
sx q[2];
rz(-2.4028845) q[2];
sx q[2];
rz(0.040741097) q[2];
rz(0.1013969) q[3];
sx q[3];
rz(-1.3359759) q[3];
sx q[3];
rz(-2.0749157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7539702) q[0];
sx q[0];
rz(-3.1356223) q[0];
sx q[0];
rz(0.23400865) q[0];
rz(-2.9435844) q[1];
sx q[1];
rz(-2.104069) q[1];
sx q[1];
rz(-2.1580946) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7957243) q[0];
sx q[0];
rz(-2.0540753) q[0];
sx q[0];
rz(-0.048075284) q[0];
x q[1];
rz(1.6621236) q[2];
sx q[2];
rz(-1.0589561) q[2];
sx q[2];
rz(2.3607766) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2605839) q[1];
sx q[1];
rz(-2.1420292) q[1];
sx q[1];
rz(-2.2991979) q[1];
rz(-pi) q[2];
rz(1.6487021) q[3];
sx q[3];
rz(-2.7512105) q[3];
sx q[3];
rz(-1.4161033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6733751) q[2];
sx q[2];
rz(-2.8523291) q[2];
sx q[2];
rz(0.77486983) q[2];
rz(-1.8153927) q[3];
sx q[3];
rz(-1.1893136) q[3];
sx q[3];
rz(2.6302122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4020017) q[0];
sx q[0];
rz(-2.2690161) q[0];
sx q[0];
rz(-0.032489754) q[0];
rz(1.2292817) q[1];
sx q[1];
rz(-2.2132497) q[1];
sx q[1];
rz(0.083018735) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5755324) q[0];
sx q[0];
rz(-1.119918) q[0];
sx q[0];
rz(1.8203446) q[0];
rz(-pi) q[1];
rz(0.55687161) q[2];
sx q[2];
rz(-0.31302127) q[2];
sx q[2];
rz(-1.0525296) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28369812) q[1];
sx q[1];
rz(-2.2304728) q[1];
sx q[1];
rz(-2.1723351) q[1];
x q[2];
rz(-3.1260499) q[3];
sx q[3];
rz(-2.2076026) q[3];
sx q[3];
rz(-1.6890845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7897537) q[2];
sx q[2];
rz(-2.3408076) q[2];
sx q[2];
rz(2.9808673) q[2];
rz(1.9488581) q[3];
sx q[3];
rz(-0.21218097) q[3];
sx q[3];
rz(-0.76782697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6650894) q[0];
sx q[0];
rz(-1.1192717) q[0];
sx q[0];
rz(-2.8506668) q[0];
rz(-1.3833969) q[1];
sx q[1];
rz(-1.8427126) q[1];
sx q[1];
rz(-0.49984041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0394589) q[0];
sx q[0];
rz(-1.5284662) q[0];
sx q[0];
rz(3.0979587) q[0];
rz(-0.96237125) q[2];
sx q[2];
rz(-0.82686868) q[2];
sx q[2];
rz(-2.8679071) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.11092) q[1];
sx q[1];
rz(-0.77092147) q[1];
sx q[1];
rz(0.67071192) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.860943) q[3];
sx q[3];
rz(-1.9925041) q[3];
sx q[3];
rz(0.56906869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9342039) q[2];
sx q[2];
rz(-0.61177212) q[2];
sx q[2];
rz(2.440051) q[2];
rz(-2.1675341) q[3];
sx q[3];
rz(-2.3249224) q[3];
sx q[3];
rz(-2.2402703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3153673) q[0];
sx q[0];
rz(-0.81804818) q[0];
sx q[0];
rz(-0.65959626) q[0];
rz(1.9024128) q[1];
sx q[1];
rz(-2.0678803) q[1];
sx q[1];
rz(-2.0674131) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91759085) q[0];
sx q[0];
rz(-1.5573643) q[0];
sx q[0];
rz(-1.55634) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42785449) q[2];
sx q[2];
rz(-1.9249467) q[2];
sx q[2];
rz(-2.4216975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2905137) q[1];
sx q[1];
rz(-2.004262) q[1];
sx q[1];
rz(0.23833935) q[1];
rz(2.6857008) q[3];
sx q[3];
rz(-2.3237202) q[3];
sx q[3];
rz(-0.65480937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2913975) q[2];
sx q[2];
rz(-0.75866282) q[2];
sx q[2];
rz(-0.58471739) q[2];
rz(-1.0662474) q[3];
sx q[3];
rz(-1.2277579) q[3];
sx q[3];
rz(-2.7339981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19121118) q[0];
sx q[0];
rz(-1.9712912) q[0];
sx q[0];
rz(2.6608652) q[0];
rz(1.9302543) q[1];
sx q[1];
rz(-1.6119266) q[1];
sx q[1];
rz(-2.5239677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6473501) q[0];
sx q[0];
rz(-1.8959008) q[0];
sx q[0];
rz(0.0044429739) q[0];
rz(2.4489538) q[2];
sx q[2];
rz(-2.4071781) q[2];
sx q[2];
rz(1.278233) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.86194) q[1];
sx q[1];
rz(-1.8649001) q[1];
sx q[1];
rz(-2.9169464) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9493773) q[3];
sx q[3];
rz(-1.6055067) q[3];
sx q[3];
rz(-0.46958967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9628613) q[2];
sx q[2];
rz(-1.3946673) q[2];
sx q[2];
rz(-2.2507131) q[2];
rz(-2.1122872) q[3];
sx q[3];
rz(-1.7512713) q[3];
sx q[3];
rz(-1.5911969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35128281) q[0];
sx q[0];
rz(-0.93515486) q[0];
sx q[0];
rz(1.9679605) q[0];
rz(1.1890746) q[1];
sx q[1];
rz(-2.3937841) q[1];
sx q[1];
rz(-0.48114166) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089398459) q[0];
sx q[0];
rz(-1.7939095) q[0];
sx q[0];
rz(0.57017489) q[0];
rz(-pi) q[1];
rz(-2.9815707) q[2];
sx q[2];
rz(-0.56567398) q[2];
sx q[2];
rz(2.6765347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4906686) q[1];
sx q[1];
rz(-1.7177525) q[1];
sx q[1];
rz(-2.1747194) q[1];
x q[2];
rz(-2.6974996) q[3];
sx q[3];
rz(-1.0837348) q[3];
sx q[3];
rz(-0.57898486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9159307) q[2];
sx q[2];
rz(-1.2292726) q[2];
sx q[2];
rz(1.4813598) q[2];
rz(0.046317421) q[3];
sx q[3];
rz(-1.1888844) q[3];
sx q[3];
rz(1.2272629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7939821) q[0];
sx q[0];
rz(-1.0337669) q[0];
sx q[0];
rz(3.0718497) q[0];
rz(-2.8522885) q[1];
sx q[1];
rz(-1.8569088) q[1];
sx q[1];
rz(2.0893673) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7523868) q[0];
sx q[0];
rz(-1.2121887) q[0];
sx q[0];
rz(1.9176656) q[0];
rz(-0.90862008) q[2];
sx q[2];
rz(-1.7025456) q[2];
sx q[2];
rz(2.3275304) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1919162) q[1];
sx q[1];
rz(-1.9468309) q[1];
sx q[1];
rz(2.7533349) q[1];
rz(-1.9931562) q[3];
sx q[3];
rz(-1.1760654) q[3];
sx q[3];
rz(2.0572061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59794402) q[2];
sx q[2];
rz(-1.1755377) q[2];
sx q[2];
rz(0.0607461) q[2];
rz(-0.28901035) q[3];
sx q[3];
rz(-1.776639) q[3];
sx q[3];
rz(1.2750767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.65171) q[0];
sx q[0];
rz(-1.4986421) q[0];
sx q[0];
rz(-1.0810252) q[0];
rz(2.091485) q[1];
sx q[1];
rz(-2.0790015) q[1];
sx q[1];
rz(1.9008295) q[1];
rz(1.8619887) q[2];
sx q[2];
rz(-2.6544901) q[2];
sx q[2];
rz(2.0915379) q[2];
rz(1.1393094) q[3];
sx q[3];
rz(-1.8284952) q[3];
sx q[3];
rz(-1.8313051) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
