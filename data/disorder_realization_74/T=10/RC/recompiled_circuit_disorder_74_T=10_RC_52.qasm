OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.4607853) q[0];
sx q[0];
rz(-2.1587125) q[0];
sx q[0];
rz(2.13184) q[0];
rz(2.9653964) q[1];
sx q[1];
rz(5.3806452) q[1];
sx q[1];
rz(7.5723958) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56611094) q[0];
sx q[0];
rz(-2.5876343) q[0];
sx q[0];
rz(-1.0004811) q[0];
rz(-pi) q[1];
rz(-2.870954) q[2];
sx q[2];
rz(-3.116313) q[2];
sx q[2];
rz(-1.6563005) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2002441) q[1];
sx q[1];
rz(-1.3625047) q[1];
sx q[1];
rz(1.8336481) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3207924) q[3];
sx q[3];
rz(-2.8651926) q[3];
sx q[3];
rz(-1.9763415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5499128) q[2];
sx q[2];
rz(-0.16273558) q[2];
sx q[2];
rz(0.13554779) q[2];
rz(0.1401976) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-0.31301096) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267589) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(-2.8813598) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(-1.3134726) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1056054) q[0];
sx q[0];
rz(-0.30788883) q[0];
sx q[0];
rz(-0.76460989) q[0];
x q[1];
rz(-2.1140695) q[2];
sx q[2];
rz(-0.55210219) q[2];
sx q[2];
rz(-3.0039624) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7105512) q[1];
sx q[1];
rz(-1.4841054) q[1];
sx q[1];
rz(-0.1608506) q[1];
rz(1.6766657) q[3];
sx q[3];
rz(-1.6459811) q[3];
sx q[3];
rz(0.74010805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.446622) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(-0.95412811) q[2];
rz(1.4387087) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057782877) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(0.69586786) q[0];
rz(2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(-0.16608873) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3375219) q[0];
sx q[0];
rz(-1.4253758) q[0];
sx q[0];
rz(1.279633) q[0];
rz(-2.0735998) q[2];
sx q[2];
rz(-1.0223801) q[2];
sx q[2];
rz(-1.3082248) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42576158) q[1];
sx q[1];
rz(-2.466723) q[1];
sx q[1];
rz(-0.75158822) q[1];
x q[2];
rz(-0.8576287) q[3];
sx q[3];
rz(-1.7664675) q[3];
sx q[3];
rz(1.5426202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7819536) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(0.98177838) q[2];
rz(-2.0180457) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(-1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(-3.0825901) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(-2.847805) q[0];
rz(-0.44149533) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-0.77484432) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9967277) q[0];
sx q[0];
rz(-1.5872247) q[0];
sx q[0];
rz(-0.073547151) q[0];
rz(-1.5261126) q[2];
sx q[2];
rz(-1.2031021) q[2];
sx q[2];
rz(-1.1553264) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6619819) q[1];
sx q[1];
rz(-1.4927215) q[1];
sx q[1];
rz(1.1675203) q[1];
rz(-2.1500548) q[3];
sx q[3];
rz(-0.30617985) q[3];
sx q[3];
rz(-1.7508208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.0094770771) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(0.57007989) q[2];
rz(-1.8360957) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(-1.501804) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.775979) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(2.9602125) q[0];
rz(0.60896215) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(3.0338874) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5126702) q[0];
sx q[0];
rz(-1.9229917) q[0];
sx q[0];
rz(-3.0432426) q[0];
rz(-pi) q[1];
rz(2.2568251) q[2];
sx q[2];
rz(-2.3713285) q[2];
sx q[2];
rz(1.2203072) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2376033) q[1];
sx q[1];
rz(-0.9346107) q[1];
sx q[1];
rz(0.49828766) q[1];
x q[2];
rz(-1.0072458) q[3];
sx q[3];
rz(-1.8661024) q[3];
sx q[3];
rz(-1.9293279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.084215) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(2.8520612) q[2];
rz(-0.036751898) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(3.0832131) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(-1.6311197) q[0];
rz(1.1389114) q[1];
sx q[1];
rz(-1.6612256) q[1];
sx q[1];
rz(2.1246134) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4879918) q[0];
sx q[0];
rz(-1.2383608) q[0];
sx q[0];
rz(1.0792653) q[0];
x q[1];
rz(-0.20118841) q[2];
sx q[2];
rz(-1.5941465) q[2];
sx q[2];
rz(0.77229283) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14156995) q[1];
sx q[1];
rz(-2.0380028) q[1];
sx q[1];
rz(1.0228553) q[1];
x q[2];
rz(-1.059504) q[3];
sx q[3];
rz(-2.0492616) q[3];
sx q[3];
rz(-2.2389776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5974474) q[2];
sx q[2];
rz(-1.0015254) q[2];
sx q[2];
rz(0.63465676) q[2];
rz(-2.0641816) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8909661) q[0];
sx q[0];
rz(-2.4996596) q[0];
sx q[0];
rz(2.2977258) q[0];
rz(-1.5232874) q[1];
sx q[1];
rz(-0.73262501) q[1];
sx q[1];
rz(1.9218146) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6780753) q[0];
sx q[0];
rz(-0.099637195) q[0];
sx q[0];
rz(-2.7479322) q[0];
x q[1];
rz(-2.8104066) q[2];
sx q[2];
rz(-1.965431) q[2];
sx q[2];
rz(-1.0489724) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0691094) q[1];
sx q[1];
rz(-1.8169889) q[1];
sx q[1];
rz(-2.2407131) q[1];
rz(-pi) q[2];
rz(3.0516902) q[3];
sx q[3];
rz(-1.9780759) q[3];
sx q[3];
rz(-2.390993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5119778) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(2.8249595) q[2];
rz(3.1043502) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(-0.75604701) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1036296) q[0];
sx q[0];
rz(-2.0541971) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(-1.9695075) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(0.40922871) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6054571) q[0];
sx q[0];
rz(-0.94263173) q[0];
sx q[0];
rz(2.7917042) q[0];
rz(2.6312469) q[2];
sx q[2];
rz(-0.58594698) q[2];
sx q[2];
rz(-2.3454391) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.65539) q[1];
sx q[1];
rz(-1.4185925) q[1];
sx q[1];
rz(0.57032077) q[1];
rz(-2.7301634) q[3];
sx q[3];
rz(-0.55555389) q[3];
sx q[3];
rz(-0.065447741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1428712) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(0.55465737) q[2];
rz(2.5000642) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91759578) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(3.1316485) q[0];
rz(-1.0558646) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(1.508629) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89440896) q[0];
sx q[0];
rz(-1.4100473) q[0];
sx q[0];
rz(-0.87662351) q[0];
rz(-pi) q[1];
rz(-3.0502003) q[2];
sx q[2];
rz(-0.95801991) q[2];
sx q[2];
rz(0.44832715) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.75084) q[1];
sx q[1];
rz(-2.5166593) q[1];
sx q[1];
rz(-1.0052488) q[1];
rz(0.86767254) q[3];
sx q[3];
rz(-1.613986) q[3];
sx q[3];
rz(1.6004576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.14279723) q[2];
sx q[2];
rz(-1.5974416) q[2];
sx q[2];
rz(-2.7049086) q[2];
rz(1.8113332) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0325539) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(0.56525266) q[0];
rz(0.25282192) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(1.7601097) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5436343) q[0];
sx q[0];
rz(-1.6126452) q[0];
sx q[0];
rz(2.6393487) q[0];
rz(-2.6673615) q[2];
sx q[2];
rz(-1.6696764) q[2];
sx q[2];
rz(2.729051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6773721) q[1];
sx q[1];
rz(-2.3090625) q[1];
sx q[1];
rz(-2.1470451) q[1];
rz(-pi) q[2];
rz(-1.2425209) q[3];
sx q[3];
rz(-2.4007113) q[3];
sx q[3];
rz(1.5433943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4505724) q[2];
sx q[2];
rz(-2.317163) q[2];
sx q[2];
rz(-2.5449469) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(3.1141282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(0.52994603) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(1.272841) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(-3.0957072) q[2];
sx q[2];
rz(-1.3165717) q[2];
sx q[2];
rz(-1.7236621) q[2];
rz(-1.6521372) q[3];
sx q[3];
rz(-1.0621536) q[3];
sx q[3];
rz(-2.5928706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
