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
rz(-1.9189605) q[0];
sx q[0];
rz(2.5834592) q[0];
sx q[0];
rz(8.7659788) q[0];
rz(0.88762033) q[1];
sx q[1];
rz(-0.69488156) q[1];
sx q[1];
rz(-3.053009) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9889116) q[0];
sx q[0];
rz(-1.7636239) q[0];
sx q[0];
rz(2.0047943) q[0];
rz(-2.1363356) q[2];
sx q[2];
rz(-1.9449678) q[2];
sx q[2];
rz(-2.1924874) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3693847) q[1];
sx q[1];
rz(-1.7225035) q[1];
sx q[1];
rz(2.1939192) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5578868) q[3];
sx q[3];
rz(-0.87432623) q[3];
sx q[3];
rz(-0.59244746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6170071) q[2];
sx q[2];
rz(-1.969939) q[2];
sx q[2];
rz(0.48801547) q[2];
rz(0.46767849) q[3];
sx q[3];
rz(-0.52507639) q[3];
sx q[3];
rz(0.16744965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92983285) q[0];
sx q[0];
rz(-2.3303895) q[0];
sx q[0];
rz(-1.2540586) q[0];
rz(-2.5414741) q[1];
sx q[1];
rz(-1.3328726) q[1];
sx q[1];
rz(-0.41807237) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9771043) q[0];
sx q[0];
rz(-0.56796861) q[0];
sx q[0];
rz(-2.8051359) q[0];
x q[1];
rz(-1.4457147) q[2];
sx q[2];
rz(-0.92601771) q[2];
sx q[2];
rz(-3.038542) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.969287) q[1];
sx q[1];
rz(-1.8864156) q[1];
sx q[1];
rz(-3.0728042) q[1];
rz(-pi) q[2];
rz(2.5097519) q[3];
sx q[3];
rz(-2.2152293) q[3];
sx q[3];
rz(2.0845678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7347001) q[2];
sx q[2];
rz(-2.0105346) q[2];
sx q[2];
rz(-0.1184173) q[2];
rz(-0.40220574) q[3];
sx q[3];
rz(-1.4879358) q[3];
sx q[3];
rz(1.3291298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6953832) q[0];
sx q[0];
rz(-0.24676794) q[0];
sx q[0];
rz(1.9870019) q[0];
rz(1.062695) q[1];
sx q[1];
rz(-1.0239536) q[1];
sx q[1];
rz(-1.9564995) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20511928) q[0];
sx q[0];
rz(-1.8471878) q[0];
sx q[0];
rz(0.87262459) q[0];
rz(-1.0555004) q[2];
sx q[2];
rz(-1.3772817) q[2];
sx q[2];
rz(-3.0119579) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1968763) q[1];
sx q[1];
rz(-1.8507974) q[1];
sx q[1];
rz(1.9619322) q[1];
rz(-pi) q[2];
rz(-0.40374741) q[3];
sx q[3];
rz(-1.8093411) q[3];
sx q[3];
rz(-0.26097894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2350754) q[2];
sx q[2];
rz(-0.86248988) q[2];
sx q[2];
rz(-1.4441351) q[2];
rz(0.92153543) q[3];
sx q[3];
rz(-0.60465616) q[3];
sx q[3];
rz(-3.0864033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.4398572) q[0];
sx q[0];
rz(-0.12049645) q[0];
sx q[0];
rz(3.0635656) q[0];
rz(2.0241418) q[1];
sx q[1];
rz(-1.2470587) q[1];
sx q[1];
rz(-1.6695581) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0271281) q[0];
sx q[0];
rz(-0.77550602) q[0];
sx q[0];
rz(-1.7239611) q[0];
x q[1];
rz(0.53773625) q[2];
sx q[2];
rz(-2.0860279) q[2];
sx q[2];
rz(-2.4071884) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8406183) q[1];
sx q[1];
rz(-1.930791) q[1];
sx q[1];
rz(3.1112413) q[1];
rz(-0.41507657) q[3];
sx q[3];
rz(-0.76722324) q[3];
sx q[3];
rz(-1.5515212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3052519) q[2];
sx q[2];
rz(-2.4938816) q[2];
sx q[2];
rz(-0.15023896) q[2];
rz(-2.5495106) q[3];
sx q[3];
rz(-1.3034857) q[3];
sx q[3];
rz(-1.9534115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93064654) q[0];
sx q[0];
rz(-0.35363126) q[0];
sx q[0];
rz(0.76552248) q[0];
rz(-2.3402479) q[1];
sx q[1];
rz(-2.0739906) q[1];
sx q[1];
rz(1.9498922) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0203637) q[0];
sx q[0];
rz(-1.0397216) q[0];
sx q[0];
rz(1.9241821) q[0];
rz(-0.72010626) q[2];
sx q[2];
rz(-2.4821447) q[2];
sx q[2];
rz(-2.5044005) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5728103) q[1];
sx q[1];
rz(-0.88925075) q[1];
sx q[1];
rz(-3.0589025) q[1];
rz(-2.3611154) q[3];
sx q[3];
rz(-2.7679518) q[3];
sx q[3];
rz(2.3821156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2661065) q[2];
sx q[2];
rz(-2.2849639) q[2];
sx q[2];
rz(0.53492707) q[2];
rz(2.5180425) q[3];
sx q[3];
rz(-2.0303191) q[3];
sx q[3];
rz(-2.258544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18735886) q[0];
sx q[0];
rz(-2.7281902) q[0];
sx q[0];
rz(-0.9340539) q[0];
rz(1.8961228) q[1];
sx q[1];
rz(-1.586069) q[1];
sx q[1];
rz(-0.62560558) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0348822) q[0];
sx q[0];
rz(-0.83256522) q[0];
sx q[0];
rz(0.24068479) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0482381) q[2];
sx q[2];
rz(-2.2389004) q[2];
sx q[2];
rz(2.7808288) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6954591) q[1];
sx q[1];
rz(-1.8605821) q[1];
sx q[1];
rz(-2.4373846) q[1];
rz(-1.7966426) q[3];
sx q[3];
rz(-1.3991688) q[3];
sx q[3];
rz(0.93935475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0921359) q[2];
sx q[2];
rz(-0.9149887) q[2];
sx q[2];
rz(0.12506872) q[2];
rz(-2.4140029) q[3];
sx q[3];
rz(-1.5743419) q[3];
sx q[3];
rz(2.6653813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26559386) q[0];
sx q[0];
rz(-0.47530526) q[0];
sx q[0];
rz(-1.884961) q[0];
rz(1.9427293) q[1];
sx q[1];
rz(-1.7190944) q[1];
sx q[1];
rz(-1.8667603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0860109) q[0];
sx q[0];
rz(-2.0025829) q[0];
sx q[0];
rz(0.55786572) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0339678) q[2];
sx q[2];
rz(-1.2481999) q[2];
sx q[2];
rz(-0.19746298) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5526455) q[1];
sx q[1];
rz(-1.1324404) q[1];
sx q[1];
rz(0.38927765) q[1];
rz(0.12673283) q[3];
sx q[3];
rz(-1.885948) q[3];
sx q[3];
rz(-2.2268068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4510497) q[2];
sx q[2];
rz(-2.1284911) q[2];
sx q[2];
rz(0.96735442) q[2];
rz(1.5585772) q[3];
sx q[3];
rz(-1.6396921) q[3];
sx q[3];
rz(-0.30174747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3485182) q[0];
sx q[0];
rz(-2.5418042) q[0];
sx q[0];
rz(0.10979688) q[0];
rz(1.9684567) q[1];
sx q[1];
rz(-0.99183142) q[1];
sx q[1];
rz(-2.0593624) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17324461) q[0];
sx q[0];
rz(-0.81615198) q[0];
sx q[0];
rz(-1.9108755) q[0];
rz(-3.033048) q[2];
sx q[2];
rz(-1.0883779) q[2];
sx q[2];
rz(-0.71297836) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.70955706) q[1];
sx q[1];
rz(-2.2974778) q[1];
sx q[1];
rz(-2.130351) q[1];
rz(1.1940597) q[3];
sx q[3];
rz(-0.96661813) q[3];
sx q[3];
rz(2.8967711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.75862306) q[2];
sx q[2];
rz(-1.0926282) q[2];
sx q[2];
rz(0.36337241) q[2];
rz(2.162497) q[3];
sx q[3];
rz(-2.4978814) q[3];
sx q[3];
rz(0.50160971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8363504) q[0];
sx q[0];
rz(-2.3452121) q[0];
sx q[0];
rz(-0.35261944) q[0];
rz(-1.7800219) q[1];
sx q[1];
rz(-1.9510599) q[1];
sx q[1];
rz(3.000066) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82194009) q[0];
sx q[0];
rz(-1.7403649) q[0];
sx q[0];
rz(-1.5362306) q[0];
rz(-2.0317215) q[2];
sx q[2];
rz(-1.0264947) q[2];
sx q[2];
rz(1.6044631) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.80672164) q[1];
sx q[1];
rz(-1.1523243) q[1];
sx q[1];
rz(1.970402) q[1];
x q[2];
rz(1.0120506) q[3];
sx q[3];
rz(-1.0025452) q[3];
sx q[3];
rz(-2.760134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3129468) q[2];
sx q[2];
rz(-1.4742278) q[2];
sx q[2];
rz(2.055577) q[2];
rz(0.18276754) q[3];
sx q[3];
rz(-1.026231) q[3];
sx q[3];
rz(0.97203794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6163841) q[0];
sx q[0];
rz(-1.5807736) q[0];
sx q[0];
rz(2.4973448) q[0];
rz(0.066990189) q[1];
sx q[1];
rz(-1.3863775) q[1];
sx q[1];
rz(3.0279874) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8380801) q[0];
sx q[0];
rz(-1.7948196) q[0];
sx q[0];
rz(0.0043078686) q[0];
rz(-pi) q[1];
rz(-0.89130104) q[2];
sx q[2];
rz(-3.0105053) q[2];
sx q[2];
rz(-2.2830414) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5731331) q[1];
sx q[1];
rz(-1.4541469) q[1];
sx q[1];
rz(-0.17165143) q[1];
rz(1.6842277) q[3];
sx q[3];
rz(-0.3598752) q[3];
sx q[3];
rz(0.42478334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99864787) q[2];
sx q[2];
rz(-1.8871658) q[2];
sx q[2];
rz(-0.76510731) q[2];
rz(0.1782002) q[3];
sx q[3];
rz(-1.5604115) q[3];
sx q[3];
rz(-2.3044738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9255623) q[0];
sx q[0];
rz(-0.75556527) q[0];
sx q[0];
rz(2.8788993) q[0];
rz(-0.33949159) q[1];
sx q[1];
rz(-1.2295634) q[1];
sx q[1];
rz(-2.0801574) q[1];
rz(0.4364261) q[2];
sx q[2];
rz(-0.70311762) q[2];
sx q[2];
rz(1.22618) q[2];
rz(-2.3332023) q[3];
sx q[3];
rz(-2.7182073) q[3];
sx q[3];
rz(-3.0207664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
