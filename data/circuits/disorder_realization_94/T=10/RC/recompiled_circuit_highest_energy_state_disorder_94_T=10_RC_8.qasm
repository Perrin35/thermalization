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
rz(1.2430159) q[0];
sx q[0];
rz(-1.1216811) q[0];
sx q[0];
rz(-2.6774874) q[0];
rz(2.3857181) q[1];
sx q[1];
rz(-0.90944374) q[1];
sx q[1];
rz(-1.3165201) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6911294) q[0];
sx q[0];
rz(-3.0678406) q[0];
sx q[0];
rz(-2.6534257) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8009802) q[2];
sx q[2];
rz(-1.5235708) q[2];
sx q[2];
rz(1.7318341) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94095549) q[1];
sx q[1];
rz(-0.66087729) q[1];
sx q[1];
rz(-2.6913068) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2795696) q[3];
sx q[3];
rz(-1.3253477) q[3];
sx q[3];
rz(-0.042118532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89062771) q[2];
sx q[2];
rz(-2.784412) q[2];
sx q[2];
rz(-2.8446021) q[2];
rz(2.4885528) q[3];
sx q[3];
rz(-1.5584471) q[3];
sx q[3];
rz(2.0853341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.5121269) q[0];
sx q[0];
rz(-1.0221721) q[0];
sx q[0];
rz(2.8566991) q[0];
rz(3.0902872) q[1];
sx q[1];
rz(-0.76301328) q[1];
sx q[1];
rz(-2.5047393) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92572407) q[0];
sx q[0];
rz(-1.1236188) q[0];
sx q[0];
rz(0.53430064) q[0];
rz(-pi) q[1];
rz(0.35164386) q[2];
sx q[2];
rz(-1.571901) q[2];
sx q[2];
rz(1.0467093) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7015345) q[1];
sx q[1];
rz(-1.4949168) q[1];
sx q[1];
rz(-0.17497356) q[1];
x q[2];
rz(1.7039744) q[3];
sx q[3];
rz(-2.2723327) q[3];
sx q[3];
rz(-3.0210036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.18179831) q[2];
sx q[2];
rz(-0.85796285) q[2];
sx q[2];
rz(-0.96724969) q[2];
rz(-2.8773384) q[3];
sx q[3];
rz(-1.5033009) q[3];
sx q[3];
rz(1.2359469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018175) q[0];
sx q[0];
rz(-1.5793261) q[0];
sx q[0];
rz(-1.1199957) q[0];
rz(-2.6657875) q[1];
sx q[1];
rz(-1.0775403) q[1];
sx q[1];
rz(0.30398223) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7426152) q[0];
sx q[0];
rz(-2.1362059) q[0];
sx q[0];
rz(-3.1088367) q[0];
x q[1];
rz(1.4690222) q[2];
sx q[2];
rz(-1.8074805) q[2];
sx q[2];
rz(-2.9477811) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7870362) q[1];
sx q[1];
rz(-0.72595412) q[1];
sx q[1];
rz(1.4158842) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7695565) q[3];
sx q[3];
rz(-2.244368) q[3];
sx q[3];
rz(-2.4186132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9756644) q[2];
sx q[2];
rz(-1.0378446) q[2];
sx q[2];
rz(-2.4456639) q[2];
rz(-1.1337229) q[3];
sx q[3];
rz(-1.1997831) q[3];
sx q[3];
rz(-0.55293647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7021779) q[0];
sx q[0];
rz(-1.004383) q[0];
sx q[0];
rz(1.0473921) q[0];
rz(1.1737431) q[1];
sx q[1];
rz(-0.67432299) q[1];
sx q[1];
rz(-1.5210927) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0899309) q[0];
sx q[0];
rz(-2.0712081) q[0];
sx q[0];
rz(-1.7080281) q[0];
rz(-pi) q[1];
rz(0.32033605) q[2];
sx q[2];
rz(-0.89010677) q[2];
sx q[2];
rz(-0.57094157) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4075267) q[1];
sx q[1];
rz(-0.95469771) q[1];
sx q[1];
rz(-2.1156838) q[1];
rz(-pi) q[2];
rz(1.8549881) q[3];
sx q[3];
rz(-1.8833369) q[3];
sx q[3];
rz(0.72060637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.095470458) q[2];
sx q[2];
rz(-2.0746524) q[2];
sx q[2];
rz(0.53786892) q[2];
rz(1.4872023) q[3];
sx q[3];
rz(-0.40707773) q[3];
sx q[3];
rz(-2.4552104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9293514) q[0];
sx q[0];
rz(-2.90726) q[0];
sx q[0];
rz(-2.5382163) q[0];
rz(-2.3790908) q[1];
sx q[1];
rz(-1.1926032) q[1];
sx q[1];
rz(0.71294436) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.514587) q[0];
sx q[0];
rz(-2.7381185) q[0];
sx q[0];
rz(1.7395942) q[0];
rz(-pi) q[1];
x q[1];
rz(1.845417) q[2];
sx q[2];
rz(-2.6634187) q[2];
sx q[2];
rz(-2.630065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.566943) q[1];
sx q[1];
rz(-1.7716367) q[1];
sx q[1];
rz(2.0748791) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6176881) q[3];
sx q[3];
rz(-1.0782982) q[3];
sx q[3];
rz(2.4503051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.9731628) q[2];
sx q[2];
rz(-0.91732401) q[2];
sx q[2];
rz(1.6707576) q[2];
rz(-2.6245608) q[3];
sx q[3];
rz(-1.0443338) q[3];
sx q[3];
rz(-1.8487336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9789199) q[0];
sx q[0];
rz(-2.6277268) q[0];
sx q[0];
rz(-2.5901219) q[0];
rz(0.035471352) q[1];
sx q[1];
rz(-1.1330117) q[1];
sx q[1];
rz(-1.5303401) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0126969) q[0];
sx q[0];
rz(-0.91980359) q[0];
sx q[0];
rz(0.76437073) q[0];
rz(-pi) q[1];
rz(0.050889767) q[2];
sx q[2];
rz(-0.82422232) q[2];
sx q[2];
rz(2.3498866) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4582461) q[1];
sx q[1];
rz(-0.85696942) q[1];
sx q[1];
rz(0.48269646) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.084613581) q[3];
sx q[3];
rz(-2.222098) q[3];
sx q[3];
rz(-2.8433329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0643206) q[2];
sx q[2];
rz(-0.52017009) q[2];
sx q[2];
rz(-1.2426097) q[2];
rz(-0.51314917) q[3];
sx q[3];
rz(-2.7276701) q[3];
sx q[3];
rz(1.7665524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4921017) q[0];
sx q[0];
rz(-0.92925564) q[0];
sx q[0];
rz(-3.0249366) q[0];
rz(0.77313441) q[1];
sx q[1];
rz(-1.9832858) q[1];
sx q[1];
rz(1.6366417) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0610675) q[0];
sx q[0];
rz(-1.8634029) q[0];
sx q[0];
rz(2.3181897) q[0];
rz(-pi) q[1];
rz(-1.9595615) q[2];
sx q[2];
rz(-1.5156399) q[2];
sx q[2];
rz(-2.9832632) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.22681397) q[1];
sx q[1];
rz(-0.89229167) q[1];
sx q[1];
rz(-2.6752276) q[1];
rz(-pi) q[2];
rz(-0.65251638) q[3];
sx q[3];
rz(-1.3400786) q[3];
sx q[3];
rz(1.1200865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92631212) q[2];
sx q[2];
rz(-0.24682385) q[2];
sx q[2];
rz(0.8872633) q[2];
rz(-2.4617646) q[3];
sx q[3];
rz(-0.66893783) q[3];
sx q[3];
rz(0.63290709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55650869) q[0];
sx q[0];
rz(-1.2932788) q[0];
sx q[0];
rz(2.4853117) q[0];
rz(0.20690021) q[1];
sx q[1];
rz(-1.2196536) q[1];
sx q[1];
rz(-1.9237178) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042699634) q[0];
sx q[0];
rz(-1.4880565) q[0];
sx q[0];
rz(-1.5356043) q[0];
rz(-pi) q[1];
rz(1.0704667) q[2];
sx q[2];
rz(-2.4917951) q[2];
sx q[2];
rz(1.4184679) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9838502) q[1];
sx q[1];
rz(-1.7291673) q[1];
sx q[1];
rz(1.6862129) q[1];
rz(-pi) q[2];
rz(-1.8402053) q[3];
sx q[3];
rz(-0.97668934) q[3];
sx q[3];
rz(1.9850933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3394341) q[2];
sx q[2];
rz(-1.755244) q[2];
sx q[2];
rz(-0.68186861) q[2];
rz(1.9631466) q[3];
sx q[3];
rz(-2.384187) q[3];
sx q[3];
rz(1.9044378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7449529) q[0];
sx q[0];
rz(-1.5322026) q[0];
sx q[0];
rz(2.9314281) q[0];
rz(-2.919803) q[1];
sx q[1];
rz(-0.91786018) q[1];
sx q[1];
rz(-3.0992357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0641441) q[0];
sx q[0];
rz(-1.3832398) q[0];
sx q[0];
rz(-1.5483039) q[0];
x q[1];
rz(-2.508411) q[2];
sx q[2];
rz(-2.0821619) q[2];
sx q[2];
rz(-0.88054576) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4952462) q[1];
sx q[1];
rz(-2.5030285) q[1];
sx q[1];
rz(0.49444838) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70019763) q[3];
sx q[3];
rz(-2.4841016) q[3];
sx q[3];
rz(-2.0746865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.28045851) q[2];
sx q[2];
rz(-2.0529842) q[2];
sx q[2];
rz(-0.66656485) q[2];
rz(2.376453) q[3];
sx q[3];
rz(-2.9355526) q[3];
sx q[3];
rz(1.4008745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6078981) q[0];
sx q[0];
rz(-0.38132897) q[0];
sx q[0];
rz(-2.4545942) q[0];
rz(-1.2835361) q[1];
sx q[1];
rz(-1.1524408) q[1];
sx q[1];
rz(2.5680465) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6015778) q[0];
sx q[0];
rz(-0.83160831) q[0];
sx q[0];
rz(2.4180146) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0230471) q[2];
sx q[2];
rz(-2.6464823) q[2];
sx q[2];
rz(-2.4972288) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4926238) q[1];
sx q[1];
rz(-2.0855806) q[1];
sx q[1];
rz(-0.3920946) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2881094) q[3];
sx q[3];
rz(-0.9694582) q[3];
sx q[3];
rz(2.7784612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0428697) q[2];
sx q[2];
rz(-1.1770153) q[2];
sx q[2];
rz(3.0549808) q[2];
rz(0.11836554) q[3];
sx q[3];
rz(-1.5990853) q[3];
sx q[3];
rz(0.043244403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31556986) q[0];
sx q[0];
rz(-2.12349) q[0];
sx q[0];
rz(0.77793599) q[0];
rz(-2.8553873) q[1];
sx q[1];
rz(-0.83120167) q[1];
sx q[1];
rz(1.3298159) q[1];
rz(-0.61997531) q[2];
sx q[2];
rz(-2.740553) q[2];
sx q[2];
rz(0.16410826) q[2];
rz(0.56747464) q[3];
sx q[3];
rz(-2.6226061) q[3];
sx q[3];
rz(-2.0593986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
