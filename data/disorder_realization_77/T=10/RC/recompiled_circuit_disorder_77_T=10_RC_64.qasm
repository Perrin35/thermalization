OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(0.16790976) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(3.0854316) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0060318) q[0];
sx q[0];
rz(-0.52838415) q[0];
sx q[0];
rz(-2.0023268) q[0];
x q[1];
rz(-0.21284717) q[2];
sx q[2];
rz(-0.93570645) q[2];
sx q[2];
rz(-2.0095306) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1311156) q[1];
sx q[1];
rz(-1.8382204) q[1];
sx q[1];
rz(1.0130151) q[1];
x q[2];
rz(-1.6151186) q[3];
sx q[3];
rz(-1.830415) q[3];
sx q[3];
rz(1.1322024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.37796676) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(-0.4326694) q[2];
rz(1.1928605) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(-0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6137961) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(1.8288076) q[0];
rz(-2.9361172) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(-1.9899433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5713455) q[0];
sx q[0];
rz(-0.76128188) q[0];
sx q[0];
rz(2.0061357) q[0];
rz(-pi) q[1];
rz(2.5580514) q[2];
sx q[2];
rz(-1.984664) q[2];
sx q[2];
rz(1.5577424) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18560219) q[1];
sx q[1];
rz(-0.62528175) q[1];
sx q[1];
rz(-0.76693265) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50206708) q[3];
sx q[3];
rz(-1.2289398) q[3];
sx q[3];
rz(1.348192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(-0.22182626) q[2];
rz(0.37718537) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8310228) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-0.036852766) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(-3.085014) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3765592) q[0];
sx q[0];
rz(-0.74900904) q[0];
sx q[0];
rz(-2.2545933) q[0];
rz(-pi) q[1];
rz(-0.47358863) q[2];
sx q[2];
rz(-2.8550672) q[2];
sx q[2];
rz(0.63602704) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6371582) q[1];
sx q[1];
rz(-1.27099) q[1];
sx q[1];
rz(-1.3624886) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9583086) q[3];
sx q[3];
rz(-2.1790677) q[3];
sx q[3];
rz(-1.1063948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8686707) q[2];
sx q[2];
rz(-1.4977095) q[2];
sx q[2];
rz(2.2154714) q[2];
rz(0.55666322) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.2264003) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(-2.3994989) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-0.4793872) q[1];
sx q[1];
rz(0.46359584) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2984943) q[0];
sx q[0];
rz(-1.8659741) q[0];
sx q[0];
rz(-0.85423268) q[0];
rz(-pi) q[1];
rz(1.5092588) q[2];
sx q[2];
rz(-1.2386285) q[2];
sx q[2];
rz(-0.57519826) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78471781) q[1];
sx q[1];
rz(-2.3349635) q[1];
sx q[1];
rz(0.23826092) q[1];
x q[2];
rz(0.47834088) q[3];
sx q[3];
rz(-2.0739177) q[3];
sx q[3];
rz(-0.92418811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7745557) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(0.049499361) q[2];
rz(-0.1285304) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(0.29770011) q[0];
rz(-2.659335) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(-2.1972426) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1221065) q[0];
sx q[0];
rz(-1.5257611) q[0];
sx q[0];
rz(1.7800063) q[0];
rz(-pi) q[1];
rz(1.96825) q[2];
sx q[2];
rz(-0.83565088) q[2];
sx q[2];
rz(3.1189001) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9152865) q[1];
sx q[1];
rz(-1.5686791) q[1];
sx q[1];
rz(1.5096942) q[1];
rz(2.6258351) q[3];
sx q[3];
rz(-2.6532647) q[3];
sx q[3];
rz(0.26667903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.258761) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(-0.10822254) q[2];
rz(-0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(-2.8172857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43679431) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(-2.5571402) q[0];
rz(-0.8862409) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(-3.086673) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6829837) q[0];
sx q[0];
rz(-0.28124547) q[0];
sx q[0];
rz(-1.2693229) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6266277) q[2];
sx q[2];
rz(-2.8179114) q[2];
sx q[2];
rz(-1.8813546) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37327787) q[1];
sx q[1];
rz(-1.0265961) q[1];
sx q[1];
rz(-1.1276223) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1543598) q[3];
sx q[3];
rz(-1.9658486) q[3];
sx q[3];
rz(2.6810255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(1.0167271) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(1.8090766) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.465437) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(-0.28453919) q[0];
rz(-2.1971205) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(-2.231266) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7401687) q[0];
sx q[0];
rz(-1.6356042) q[0];
sx q[0];
rz(3.0868953) q[0];
rz(-pi) q[1];
rz(-1.1548642) q[2];
sx q[2];
rz(-1.8850733) q[2];
sx q[2];
rz(1.2736125) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4714204) q[1];
sx q[1];
rz(-2.6176665) q[1];
sx q[1];
rz(-2.7736204) q[1];
rz(-pi) q[2];
rz(0.72083731) q[3];
sx q[3];
rz(-2.0257054) q[3];
sx q[3];
rz(2.7618559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(-0.92203036) q[2];
rz(-0.56728029) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(-1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(0.63240504) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(-2.8410889) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7104615) q[0];
sx q[0];
rz(-2.3346402) q[0];
sx q[0];
rz(1.0459082) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95327611) q[2];
sx q[2];
rz(-0.91070181) q[2];
sx q[2];
rz(2.2613139) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1860736) q[1];
sx q[1];
rz(-1.3961853) q[1];
sx q[1];
rz(-0.34592918) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31605966) q[3];
sx q[3];
rz(-0.83871597) q[3];
sx q[3];
rz(2.4162606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.58632103) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-2.3596181) q[2];
rz(2.590495) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(-2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(-3.0138299) q[0];
rz(2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(-2.382747) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42620537) q[0];
sx q[0];
rz(-1.1466768) q[0];
sx q[0];
rz(-3.12294) q[0];
x q[1];
rz(-1.7636289) q[2];
sx q[2];
rz(-0.97059965) q[2];
sx q[2];
rz(-0.45229518) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2821741) q[1];
sx q[1];
rz(-0.71600435) q[1];
sx q[1];
rz(0.23846682) q[1];
rz(1.5522478) q[3];
sx q[3];
rz(-2.4759001) q[3];
sx q[3];
rz(-1.1902283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0163429) q[2];
sx q[2];
rz(-1.3602076) q[2];
sx q[2];
rz(-0.49003595) q[2];
rz(1.7193433) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(-1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816417) q[0];
sx q[0];
rz(-0.61976969) q[0];
sx q[0];
rz(3.066257) q[0];
rz(-0.8967337) q[1];
sx q[1];
rz(-1.9672085) q[1];
sx q[1];
rz(-0.60992253) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7216435) q[0];
sx q[0];
rz(-2.3099265) q[0];
sx q[0];
rz(1.9589817) q[0];
rz(-pi) q[1];
rz(1.8462734) q[2];
sx q[2];
rz(-1.9343978) q[2];
sx q[2];
rz(-1.7737349) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8049106) q[1];
sx q[1];
rz(-1.33178) q[1];
sx q[1];
rz(2.3906624) q[1];
x q[2];
rz(-2.6033953) q[3];
sx q[3];
rz(-2.3231069) q[3];
sx q[3];
rz(-1.324211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.909409) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(2.4278736) q[2];
rz(-0.37832007) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(-0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338035) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(-2.4907885) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(-1.3748319) q[2];
sx q[2];
rz(-1.5407731) q[2];
sx q[2];
rz(-2.4827448) q[2];
rz(2.8888632) q[3];
sx q[3];
rz(-2.0287632) q[3];
sx q[3];
rz(-0.91026929) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
