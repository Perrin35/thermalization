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
rz(2.3447073) q[0];
sx q[0];
rz(-1.7303884) q[0];
sx q[0];
rz(-0.91941961) q[0];
rz(-1.8384276) q[1];
sx q[1];
rz(-3.0923831) q[1];
sx q[1];
rz(-2.4278909) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.365276) q[0];
sx q[0];
rz(-0.31319022) q[0];
sx q[0];
rz(-1.0745144) q[0];
rz(-pi) q[1];
rz(-3.0786985) q[2];
sx q[2];
rz(-1.6990635) q[2];
sx q[2];
rz(-0.72606444) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9992396) q[1];
sx q[1];
rz(-2.0217784) q[1];
sx q[1];
rz(1.5692479) q[1];
rz(-1.6380129) q[3];
sx q[3];
rz(-0.47223642) q[3];
sx q[3];
rz(-0.41154644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.087622341) q[2];
sx q[2];
rz(-2.7555608) q[2];
sx q[2];
rz(-2.7719882) q[2];
rz(1.1450279) q[3];
sx q[3];
rz(-1.6686882) q[3];
sx q[3];
rz(-2.7692774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10577781) q[0];
sx q[0];
rz(-1.6749629) q[0];
sx q[0];
rz(2.5685487) q[0];
rz(-2.0447958) q[1];
sx q[1];
rz(-1.6005452) q[1];
sx q[1];
rz(2.6699452) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99804634) q[0];
sx q[0];
rz(-2.5295288) q[0];
sx q[0];
rz(1.3785115) q[0];
x q[1];
rz(0.69219442) q[2];
sx q[2];
rz(-2.080215) q[2];
sx q[2];
rz(-2.8329527) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0309033) q[1];
sx q[1];
rz(-0.86319268) q[1];
sx q[1];
rz(-0.75548197) q[1];
x q[2];
rz(0.92949683) q[3];
sx q[3];
rz(-2.5166582) q[3];
sx q[3];
rz(2.9004824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4196709) q[2];
sx q[2];
rz(-2.2025043) q[2];
sx q[2];
rz(-2.0527077) q[2];
rz(-1.0645083) q[3];
sx q[3];
rz(-2.0239315) q[3];
sx q[3];
rz(-0.3256807) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0786781) q[0];
sx q[0];
rz(-0.18904541) q[0];
sx q[0];
rz(-0.017070008) q[0];
rz(2.4315289) q[1];
sx q[1];
rz(-2.3080669) q[1];
sx q[1];
rz(0.56627083) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44646548) q[0];
sx q[0];
rz(-2.4003865) q[0];
sx q[0];
rz(1.5873853) q[0];
rz(-2.1115569) q[2];
sx q[2];
rz(-2.2077201) q[2];
sx q[2];
rz(-0.43535296) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6614187) q[1];
sx q[1];
rz(-1.6557449) q[1];
sx q[1];
rz(0.18710356) q[1];
rz(-0.2416824) q[3];
sx q[3];
rz(-0.2727601) q[3];
sx q[3];
rz(-0.54995104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8038586) q[2];
sx q[2];
rz(-2.2457819) q[2];
sx q[2];
rz(0.0090948661) q[2];
rz(0.13633063) q[3];
sx q[3];
rz(-2.3970042) q[3];
sx q[3];
rz(-1.2135308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71753865) q[0];
sx q[0];
rz(-2.461705) q[0];
sx q[0];
rz(0.16615443) q[0];
rz(-1.1135788) q[1];
sx q[1];
rz(-0.49003092) q[1];
sx q[1];
rz(2.9794433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1095162) q[0];
sx q[0];
rz(-1.4776395) q[0];
sx q[0];
rz(1.7514125) q[0];
x q[1];
rz(-0.040302868) q[2];
sx q[2];
rz(-2.5922814) q[2];
sx q[2];
rz(1.9343164) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6238193) q[1];
sx q[1];
rz(-0.86408593) q[1];
sx q[1];
rz(-2.049974) q[1];
rz(0.6643296) q[3];
sx q[3];
rz(-1.2264612) q[3];
sx q[3];
rz(-0.20928247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1589511) q[2];
sx q[2];
rz(-2.0621767) q[2];
sx q[2];
rz(2.6137433) q[2];
rz(0.85150254) q[3];
sx q[3];
rz(-0.32367555) q[3];
sx q[3];
rz(0.94312704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8714137) q[0];
sx q[0];
rz(-1.5988007) q[0];
sx q[0];
rz(-0.80108368) q[0];
rz(0.78394765) q[1];
sx q[1];
rz(-0.74257094) q[1];
sx q[1];
rz(-0.98091006) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50612569) q[0];
sx q[0];
rz(-2.6775595) q[0];
sx q[0];
rz(2.6543174) q[0];
rz(0.37600125) q[2];
sx q[2];
rz(-1.094295) q[2];
sx q[2];
rz(0.27793542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0438761) q[1];
sx q[1];
rz(-1.4773158) q[1];
sx q[1];
rz(3.082117) q[1];
rz(-pi) q[2];
rz(-1.1889691) q[3];
sx q[3];
rz(-0.5853879) q[3];
sx q[3];
rz(-0.25661925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.058978) q[2];
sx q[2];
rz(-1.5951944) q[2];
sx q[2];
rz(0.8832461) q[2];
rz(0.20032459) q[3];
sx q[3];
rz(-0.89503461) q[3];
sx q[3];
rz(2.0506355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14466318) q[0];
sx q[0];
rz(-1.0522333) q[0];
sx q[0];
rz(-2.1494179) q[0];
rz(-2.3400173) q[1];
sx q[1];
rz(-1.293332) q[1];
sx q[1];
rz(1.7209524) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0376037) q[0];
sx q[0];
rz(-0.9524571) q[0];
sx q[0];
rz(-1.9644587) q[0];
rz(-pi) q[1];
rz(2.1031441) q[2];
sx q[2];
rz(-1.6611929) q[2];
sx q[2];
rz(0.20648512) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3657377) q[1];
sx q[1];
rz(-0.72682646) q[1];
sx q[1];
rz(-2.6353542) q[1];
rz(-pi) q[2];
rz(-1.9377557) q[3];
sx q[3];
rz(-3.0375069) q[3];
sx q[3];
rz(2.6534125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8617323) q[2];
sx q[2];
rz(-1.4272775) q[2];
sx q[2];
rz(-0.53708616) q[2];
rz(-0.01072695) q[3];
sx q[3];
rz(-2.3948632) q[3];
sx q[3];
rz(0.85095325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.2456197) q[0];
sx q[0];
rz(-2.8697822) q[0];
sx q[0];
rz(-0.33988345) q[0];
rz(-1.8396359) q[1];
sx q[1];
rz(-2.1959031) q[1];
sx q[1];
rz(-1.1713015) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7599277) q[0];
sx q[0];
rz(-1.5345083) q[0];
sx q[0];
rz(-0.29706232) q[0];
rz(-pi) q[1];
rz(-1.7843855) q[2];
sx q[2];
rz(-1.3554975) q[2];
sx q[2];
rz(-0.70206308) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9885607) q[1];
sx q[1];
rz(-1.5383729) q[1];
sx q[1];
rz(0.02213636) q[1];
rz(-pi) q[2];
rz(2.9006935) q[3];
sx q[3];
rz(-2.5031075) q[3];
sx q[3];
rz(1.0257667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12477144) q[2];
sx q[2];
rz(-1.7243959) q[2];
sx q[2];
rz(0.67145124) q[2];
rz(-0.5365544) q[3];
sx q[3];
rz(-2.1814929) q[3];
sx q[3];
rz(1.051739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3143828) q[0];
sx q[0];
rz(-2.0874513) q[0];
sx q[0];
rz(2.2453454) q[0];
rz(2.8889636) q[1];
sx q[1];
rz(-0.8600421) q[1];
sx q[1];
rz(-1.0505229) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3504068) q[0];
sx q[0];
rz(-1.6210272) q[0];
sx q[0];
rz(-2.6149261) q[0];
rz(-pi) q[1];
rz(1.7017548) q[2];
sx q[2];
rz(-2.2316861) q[2];
sx q[2];
rz(0.63643989) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.72880367) q[1];
sx q[1];
rz(-2.0319684) q[1];
sx q[1];
rz(2.9739266) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0862634) q[3];
sx q[3];
rz(-0.75859514) q[3];
sx q[3];
rz(2.7531695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5184021) q[2];
sx q[2];
rz(-0.16885997) q[2];
sx q[2];
rz(1.5625578) q[2];
rz(-1.3671499) q[3];
sx q[3];
rz(-2.0953777) q[3];
sx q[3];
rz(2.3794543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6147181) q[0];
sx q[0];
rz(-1.389483) q[0];
sx q[0];
rz(-2.7516464) q[0];
rz(-0.82603106) q[1];
sx q[1];
rz(-0.69458687) q[1];
sx q[1];
rz(2.0894076) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74653502) q[0];
sx q[0];
rz(-1.2647226) q[0];
sx q[0];
rz(2.2827882) q[0];
x q[1];
rz(-2.7446943) q[2];
sx q[2];
rz(-1.7328615) q[2];
sx q[2];
rz(1.4805178) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4349367) q[1];
sx q[1];
rz(-1.2806007) q[1];
sx q[1];
rz(1.9761843) q[1];
rz(-pi) q[2];
rz(-0.94433208) q[3];
sx q[3];
rz(-1.9667224) q[3];
sx q[3];
rz(-0.79352165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0803926) q[2];
sx q[2];
rz(-1.8597417) q[2];
sx q[2];
rz(-2.6123987) q[2];
rz(0.75096327) q[3];
sx q[3];
rz(-0.96143985) q[3];
sx q[3];
rz(-0.19415893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26218364) q[0];
sx q[0];
rz(-0.59469596) q[0];
sx q[0];
rz(1.9770812) q[0];
rz(-1.2871845) q[1];
sx q[1];
rz(-2.4559805) q[1];
sx q[1];
rz(-2.6374292) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2874448) q[0];
sx q[0];
rz(-1.7896255) q[0];
sx q[0];
rz(-1.3059421) q[0];
rz(2.7648534) q[2];
sx q[2];
rz(-1.552887) q[2];
sx q[2];
rz(-2.0359382) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87889844) q[1];
sx q[1];
rz(-2.8126908) q[1];
sx q[1];
rz(-0.94543381) q[1];
rz(2.9243117) q[3];
sx q[3];
rz(-0.8978399) q[3];
sx q[3];
rz(-2.6126044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9265499) q[2];
sx q[2];
rz(-1.5861009) q[2];
sx q[2];
rz(-3.1070993) q[2];
rz(-1.5283594) q[3];
sx q[3];
rz(-0.72051636) q[3];
sx q[3];
rz(2.8300986) q[3];
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
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16019776) q[0];
sx q[0];
rz(-1.5481411) q[0];
sx q[0];
rz(2.7465469) q[0];
rz(2.1389217) q[1];
sx q[1];
rz(-1.6802588) q[1];
sx q[1];
rz(-0.37793876) q[1];
rz(2.7591697) q[2];
sx q[2];
rz(-2.6527846) q[2];
sx q[2];
rz(-2.6058886) q[2];
rz(-0.18319753) q[3];
sx q[3];
rz(-2.8092794) q[3];
sx q[3];
rz(2.487779) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
