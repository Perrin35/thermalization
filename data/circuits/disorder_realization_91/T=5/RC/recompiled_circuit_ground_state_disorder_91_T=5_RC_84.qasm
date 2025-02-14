OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.087091669) q[0];
sx q[0];
rz(-2.2890685) q[0];
sx q[0];
rz(-0.26242119) q[0];
rz(2.8334795) q[1];
sx q[1];
rz(-2.3051655) q[1];
sx q[1];
rz(-0.0094553789) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1605837) q[0];
sx q[0];
rz(-2.0045337) q[0];
sx q[0];
rz(-1.0354614) q[0];
x q[1];
rz(1.0105825) q[2];
sx q[2];
rz(-0.081801266) q[2];
sx q[2];
rz(-2.8103925) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5061599) q[1];
sx q[1];
rz(-2.3316158) q[1];
sx q[1];
rz(0.86829888) q[1];
rz(0.81266021) q[3];
sx q[3];
rz(-1.9403807) q[3];
sx q[3];
rz(-0.11358914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0414163) q[2];
sx q[2];
rz(-2.2569423) q[2];
sx q[2];
rz(0.46119383) q[2];
rz(-2.6395116) q[3];
sx q[3];
rz(-2.3177948) q[3];
sx q[3];
rz(1.4095149) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0111888) q[0];
sx q[0];
rz(-1.6809502) q[0];
sx q[0];
rz(-2.9050264) q[0];
rz(-2.323281) q[1];
sx q[1];
rz(-1.2032443) q[1];
sx q[1];
rz(0.94258211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0088975178) q[0];
sx q[0];
rz(-0.030936154) q[0];
sx q[0];
rz(1.4832458) q[0];
rz(-pi) q[1];
rz(-0.06748345) q[2];
sx q[2];
rz(-1.6915671) q[2];
sx q[2];
rz(1.6492594) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68028223) q[1];
sx q[1];
rz(-2.4066917) q[1];
sx q[1];
rz(-3.1394385) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34063142) q[3];
sx q[3];
rz(-2.0320722) q[3];
sx q[3];
rz(-3.1203062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73611034) q[2];
sx q[2];
rz(-2.5030899) q[2];
sx q[2];
rz(-2.45641) q[2];
rz(-2.3403366) q[3];
sx q[3];
rz(-1.6359436) q[3];
sx q[3];
rz(1.8906458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291229) q[0];
sx q[0];
rz(-0.47054371) q[0];
sx q[0];
rz(-2.1614918) q[0];
rz(-0.12067548) q[1];
sx q[1];
rz(-1.9735034) q[1];
sx q[1];
rz(1.4792222) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2421766) q[0];
sx q[0];
rz(-2.3606395) q[0];
sx q[0];
rz(0.23069464) q[0];
rz(-pi) q[1];
rz(-1.1086039) q[2];
sx q[2];
rz(-2.0304907) q[2];
sx q[2];
rz(-0.60079702) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.90872902) q[1];
sx q[1];
rz(-1.5655978) q[1];
sx q[1];
rz(3.1326816) q[1];
x q[2];
rz(2.6605786) q[3];
sx q[3];
rz(-2.4949007) q[3];
sx q[3];
rz(2.8493136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7145308) q[2];
sx q[2];
rz(-1.4947299) q[2];
sx q[2];
rz(-2.96116) q[2];
rz(-2.7212972) q[3];
sx q[3];
rz(-0.97729483) q[3];
sx q[3];
rz(2.596627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9412398) q[0];
sx q[0];
rz(-1.7914597) q[0];
sx q[0];
rz(-0.049276503) q[0];
rz(-0.73206466) q[1];
sx q[1];
rz(-0.8492291) q[1];
sx q[1];
rz(1.3759618) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9889209) q[0];
sx q[0];
rz(-1.7537358) q[0];
sx q[0];
rz(1.4519917) q[0];
rz(-1.3591197) q[2];
sx q[2];
rz(-0.81640252) q[2];
sx q[2];
rz(2.8704081) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0323309) q[1];
sx q[1];
rz(-1.9606661) q[1];
sx q[1];
rz(-1.6435739) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2985372) q[3];
sx q[3];
rz(-0.84765654) q[3];
sx q[3];
rz(-1.3031194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0338318) q[2];
sx q[2];
rz(-1.7511448) q[2];
sx q[2];
rz(0.66229171) q[2];
rz(1.3307339) q[3];
sx q[3];
rz(-2.3190506) q[3];
sx q[3];
rz(-1.2983769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0996284) q[0];
sx q[0];
rz(-2.8350416) q[0];
sx q[0];
rz(1.0076667) q[0];
rz(-3.0362466) q[1];
sx q[1];
rz(-2.4844929) q[1];
sx q[1];
rz(-2.9109921) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661575) q[0];
sx q[0];
rz(-0.49336067) q[0];
sx q[0];
rz(2.3301279) q[0];
rz(-0.23955524) q[2];
sx q[2];
rz(-1.2192246) q[2];
sx q[2];
rz(-2.7186175) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8457875) q[1];
sx q[1];
rz(-1.485386) q[1];
sx q[1];
rz(-2.9763362) q[1];
rz(-pi) q[2];
rz(1.4197639) q[3];
sx q[3];
rz(-2.2103035) q[3];
sx q[3];
rz(-3.1139656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0577804) q[2];
sx q[2];
rz(-0.96692204) q[2];
sx q[2];
rz(2.9528565) q[2];
rz(-1.2356637) q[3];
sx q[3];
rz(-2.6344447) q[3];
sx q[3];
rz(-1.2602826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1756303) q[0];
sx q[0];
rz(-1.9251134) q[0];
sx q[0];
rz(2.8438582) q[0];
rz(3.0689902) q[1];
sx q[1];
rz(-2.4563792) q[1];
sx q[1];
rz(-2.4593478) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52132359) q[0];
sx q[0];
rz(-0.96858874) q[0];
sx q[0];
rz(0.35767718) q[0];
rz(-0.14839006) q[2];
sx q[2];
rz(-0.94252693) q[2];
sx q[2];
rz(-0.65199696) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.161769) q[1];
sx q[1];
rz(-1.7186972) q[1];
sx q[1];
rz(-1.921341) q[1];
rz(-2.7704084) q[3];
sx q[3];
rz(-2.1371578) q[3];
sx q[3];
rz(-1.0584099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5630774) q[2];
sx q[2];
rz(-0.35334057) q[2];
sx q[2];
rz(0.087470857) q[2];
rz(0.89720094) q[3];
sx q[3];
rz(-1.8950491) q[3];
sx q[3];
rz(2.7520666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0422269) q[0];
sx q[0];
rz(-2.0392188) q[0];
sx q[0];
rz(-1.1627831) q[0];
rz(2.9258264) q[1];
sx q[1];
rz(-0.81753221) q[1];
sx q[1];
rz(-0.93528265) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.640808) q[0];
sx q[0];
rz(-1.3324059) q[0];
sx q[0];
rz(2.3604908) q[0];
rz(0.14752702) q[2];
sx q[2];
rz(-1.8263683) q[2];
sx q[2];
rz(0.35596656) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4295807) q[1];
sx q[1];
rz(-1.64974) q[1];
sx q[1];
rz(-0.78671766) q[1];
rz(-2.8264826) q[3];
sx q[3];
rz(-1.2602046) q[3];
sx q[3];
rz(2.9950708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0995348) q[2];
sx q[2];
rz(-2.4877986) q[2];
sx q[2];
rz(2.330244) q[2];
rz(-0.30284303) q[3];
sx q[3];
rz(-1.1648014) q[3];
sx q[3];
rz(2.5109049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.621156) q[0];
sx q[0];
rz(-0.32670894) q[0];
sx q[0];
rz(-0.69123554) q[0];
rz(2.4825545) q[1];
sx q[1];
rz(-1.5182511) q[1];
sx q[1];
rz(-2.5504327) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44849075) q[0];
sx q[0];
rz(-1.8536942) q[0];
sx q[0];
rz(3.1295524) q[0];
rz(-pi) q[1];
rz(-2.8560258) q[2];
sx q[2];
rz(-1.0896249) q[2];
sx q[2];
rz(-2.8701631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2745251) q[1];
sx q[1];
rz(-0.82885107) q[1];
sx q[1];
rz(-1.1222003) q[1];
rz(-3.0847315) q[3];
sx q[3];
rz(-1.2788075) q[3];
sx q[3];
rz(-2.9546471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9968694) q[2];
sx q[2];
rz(-0.99088061) q[2];
sx q[2];
rz(1.4608176) q[2];
rz(2.3878494) q[3];
sx q[3];
rz(-1.6921348) q[3];
sx q[3];
rz(-0.47372216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.6853471) q[0];
sx q[0];
rz(-2.0401177) q[0];
sx q[0];
rz(-2.9154678) q[0];
rz(1.6926758) q[1];
sx q[1];
rz(-0.96335226) q[1];
sx q[1];
rz(1.0015782) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9739363) q[0];
sx q[0];
rz(-1.4261803) q[0];
sx q[0];
rz(-0.58982106) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6993352) q[2];
sx q[2];
rz(-1.275389) q[2];
sx q[2];
rz(-2.4396074) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16538936) q[1];
sx q[1];
rz(-1.3829297) q[1];
sx q[1];
rz(-1.6783183) q[1];
x q[2];
rz(-1.8705192) q[3];
sx q[3];
rz(-2.2329438) q[3];
sx q[3];
rz(-1.7598835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61233026) q[2];
sx q[2];
rz(-2.3713106) q[2];
sx q[2];
rz(0.15753499) q[2];
rz(-0.15466386) q[3];
sx q[3];
rz(-1.1636584) q[3];
sx q[3];
rz(0.4304339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4100819) q[0];
sx q[0];
rz(-1.9022576) q[0];
sx q[0];
rz(-2.4391158) q[0];
rz(-2.6055873) q[1];
sx q[1];
rz(-0.82967007) q[1];
sx q[1];
rz(1.747267) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4634906) q[0];
sx q[0];
rz(-2.9759088) q[0];
sx q[0];
rz(-1.413762) q[0];
x q[1];
rz(-1.6800266) q[2];
sx q[2];
rz(-1.6642206) q[2];
sx q[2];
rz(-1.0051073) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.179217) q[1];
sx q[1];
rz(-1.0581279) q[1];
sx q[1];
rz(2.8151399) q[1];
x q[2];
rz(0.0004041632) q[3];
sx q[3];
rz(-0.41557103) q[3];
sx q[3];
rz(0.11041752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.006762) q[2];
sx q[2];
rz(-2.4595021) q[2];
sx q[2];
rz(-1.4466059) q[2];
rz(0.26099482) q[3];
sx q[3];
rz(-2.1311396) q[3];
sx q[3];
rz(-1.235777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21988729) q[0];
sx q[0];
rz(-2.4867262) q[0];
sx q[0];
rz(2.126271) q[0];
rz(1.075853) q[1];
sx q[1];
rz(-1.8668108) q[1];
sx q[1];
rz(1.6783953) q[1];
rz(0.27412065) q[2];
sx q[2];
rz(-0.39015301) q[2];
sx q[2];
rz(2.9377666) q[2];
rz(2.0113284) q[3];
sx q[3];
rz(-2.2838099) q[3];
sx q[3];
rz(1.6259307) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
