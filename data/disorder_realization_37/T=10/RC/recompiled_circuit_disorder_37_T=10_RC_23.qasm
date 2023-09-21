OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(-1.9394983) q[0];
sx q[0];
rz(-1.1480968) q[0];
rz(1.2530874) q[1];
sx q[1];
rz(4.0822786) q[1];
sx q[1];
rz(10.818426) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049126547) q[0];
sx q[0];
rz(-1.2611715) q[0];
sx q[0];
rz(1.5686839) q[0];
rz(-pi) q[1];
rz(0.78656466) q[2];
sx q[2];
rz(-2.2872891) q[2];
sx q[2];
rz(-0.63209817) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4268036) q[1];
sx q[1];
rz(-2.0232632) q[1];
sx q[1];
rz(1.8718375) q[1];
rz(2.1288793) q[3];
sx q[3];
rz(-1.8126243) q[3];
sx q[3];
rz(-1.1954607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.48774886) q[2];
sx q[2];
rz(-1.2922492) q[2];
sx q[2];
rz(0.1208819) q[2];
rz(-0.17928784) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.091846175) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(0.13277408) q[0];
rz(1.6800539) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(0.24138385) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5854908) q[0];
sx q[0];
rz(-1.1855159) q[0];
sx q[0];
rz(2.7566064) q[0];
rz(-pi) q[1];
rz(2.0627459) q[2];
sx q[2];
rz(-1.098512) q[2];
sx q[2];
rz(-0.3837331) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1387716) q[1];
sx q[1];
rz(-1.661794) q[1];
sx q[1];
rz(-0.9072733) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4003795) q[3];
sx q[3];
rz(-2.9885871) q[3];
sx q[3];
rz(-1.1982329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0138578) q[2];
sx q[2];
rz(-1.8227791) q[2];
sx q[2];
rz(2.0347118) q[2];
rz(1.3876623) q[3];
sx q[3];
rz(-2.6219086) q[3];
sx q[3];
rz(-1.2319516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7805507) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(-0.74044359) q[0];
rz(2.6904147) q[1];
sx q[1];
rz(-1.2606882) q[1];
sx q[1];
rz(-1.0528475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42441472) q[0];
sx q[0];
rz(-2.4632235) q[0];
sx q[0];
rz(0.60027392) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29849507) q[2];
sx q[2];
rz(-0.91909354) q[2];
sx q[2];
rz(-3.0534844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14245089) q[1];
sx q[1];
rz(-1.4186727) q[1];
sx q[1];
rz(1.4494447) q[1];
x q[2];
rz(-0.0099817688) q[3];
sx q[3];
rz(-2.0008893) q[3];
sx q[3];
rz(0.86356589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.25049245) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(0.54692522) q[2];
rz(0.28918239) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(-0.012399013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1716487) q[0];
sx q[0];
rz(-0.26370731) q[0];
sx q[0];
rz(1.3522211) q[0];
rz(2.9317454) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(2.9052177) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0683595) q[0];
sx q[0];
rz(-1.3993565) q[0];
sx q[0];
rz(1.839848) q[0];
rz(-pi) q[1];
x q[1];
rz(1.054748) q[2];
sx q[2];
rz(-1.8785254) q[2];
sx q[2];
rz(-0.2211406) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.216815) q[1];
sx q[1];
rz(-1.5642011) q[1];
sx q[1];
rz(1.4031246) q[1];
rz(-pi) q[2];
rz(1.719791) q[3];
sx q[3];
rz(-2.0553203) q[3];
sx q[3];
rz(1.5414343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9005047) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(-0.55348712) q[2];
rz(-2.2262946) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6699162) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(-2.4374403) q[0];
rz(2.0856805) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(2.5456837) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0376315) q[0];
sx q[0];
rz(-0.2346633) q[0];
sx q[0];
rz(-1.1681359) q[0];
x q[1];
rz(0.67363588) q[2];
sx q[2];
rz(-1.6792225) q[2];
sx q[2];
rz(3.0631531) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.039918598) q[1];
sx q[1];
rz(-2.6933751) q[1];
sx q[1];
rz(-0.43804534) q[1];
rz(-pi) q[2];
rz(-1.7549873) q[3];
sx q[3];
rz(-1.9362238) q[3];
sx q[3];
rz(-2.3900677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.492505) q[2];
sx q[2];
rz(-2.0272144) q[2];
sx q[2];
rz(-2.5904783) q[2];
rz(-0.20714949) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15774396) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(-2.1898848) q[0];
rz(0.47479182) q[1];
sx q[1];
rz(-1.1750849) q[1];
sx q[1];
rz(0.29528433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6280569) q[0];
sx q[0];
rz(-0.4204458) q[0];
sx q[0];
rz(-1.9168617) q[0];
x q[1];
rz(-2.6906602) q[2];
sx q[2];
rz(-0.34313289) q[2];
sx q[2];
rz(-2.1692838) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.59776781) q[1];
sx q[1];
rz(-2.4104767) q[1];
sx q[1];
rz(0.75146971) q[1];
x q[2];
rz(0.065159273) q[3];
sx q[3];
rz(-2.0805196) q[3];
sx q[3];
rz(2.5173957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39020145) q[2];
sx q[2];
rz(-1.3662806) q[2];
sx q[2];
rz(-0.93377101) q[2];
rz(-1.4592524) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(-1.6850083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0657601) q[0];
sx q[0];
rz(-1.9969143) q[0];
sx q[0];
rz(0.24205762) q[0];
rz(-0.66479713) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(0.40329969) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4664073) q[0];
sx q[0];
rz(-1.2012321) q[0];
sx q[0];
rz(0.36375605) q[0];
rz(-pi) q[1];
rz(-1.2861916) q[2];
sx q[2];
rz(-0.91214123) q[2];
sx q[2];
rz(1.8919485) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1937716) q[1];
sx q[1];
rz(-1.3606451) q[1];
sx q[1];
rz(0.0019046849) q[1];
rz(1.7607479) q[3];
sx q[3];
rz(-1.5459832) q[3];
sx q[3];
rz(-2.1638889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2287801) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(-2.7071803) q[2];
rz(-1.0007535) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(2.6312857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1338761) q[0];
sx q[0];
rz(-3.1354597) q[0];
sx q[0];
rz(0.49466053) q[0];
rz(1.5085295) q[1];
sx q[1];
rz(-1.6260908) q[1];
sx q[1];
rz(-0.60044926) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8884044) q[0];
sx q[0];
rz(-1.1115523) q[0];
sx q[0];
rz(1.123239) q[0];
rz(-pi) q[1];
rz(-1.0457439) q[2];
sx q[2];
rz(-1.5935491) q[2];
sx q[2];
rz(-2.3843228) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9143876) q[1];
sx q[1];
rz(-1.3890508) q[1];
sx q[1];
rz(1.449613) q[1];
rz(-pi) q[2];
rz(0.59372254) q[3];
sx q[3];
rz(-2.7799118) q[3];
sx q[3];
rz(0.42632521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1373458) q[2];
sx q[2];
rz(-2.1477951) q[2];
sx q[2];
rz(0.11631575) q[2];
rz(-0.4256734) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(-1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.9882934) q[0];
sx q[0];
rz(-0.17803742) q[0];
sx q[0];
rz(-1.4784038) q[0];
rz(-0.93961811) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(-0.41752648) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5719713) q[0];
sx q[0];
rz(-2.6082391) q[0];
sx q[0];
rz(1.6589952) q[0];
x q[1];
rz(-1.9475627) q[2];
sx q[2];
rz(-2.8586839) q[2];
sx q[2];
rz(0.49809581) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7716277) q[1];
sx q[1];
rz(-1.4039478) q[1];
sx q[1];
rz(-0.19373993) q[1];
rz(-1.1674676) q[3];
sx q[3];
rz(-0.46357337) q[3];
sx q[3];
rz(1.7183255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6614762) q[2];
sx q[2];
rz(-2.5051703) q[2];
sx q[2];
rz(2.647906) q[2];
rz(2.4168329) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(0.31989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9676554) q[0];
sx q[0];
rz(-0.65615654) q[0];
sx q[0];
rz(-0.68558145) q[0];
rz(-0.29742345) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(-1.1313653) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.92537) q[0];
sx q[0];
rz(-0.67514456) q[0];
sx q[0];
rz(2.569909) q[0];
x q[1];
rz(-0.71815021) q[2];
sx q[2];
rz(-2.2946723) q[2];
sx q[2];
rz(1.3494327) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.28045052) q[1];
sx q[1];
rz(-0.44499731) q[1];
sx q[1];
rz(-2.599237) q[1];
x q[2];
rz(-2.0496619) q[3];
sx q[3];
rz(-1.0189971) q[3];
sx q[3];
rz(2.5294876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.82548213) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(2.4342009) q[2];
rz(-0.80983821) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(-0.18856089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903044) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(-2.9293625) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(2.390776) q[2];
sx q[2];
rz(-2.3859947) q[2];
sx q[2];
rz(1.044556) q[2];
rz(-0.57701941) q[3];
sx q[3];
rz(-2.1885625) q[3];
sx q[3];
rz(2.267005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];