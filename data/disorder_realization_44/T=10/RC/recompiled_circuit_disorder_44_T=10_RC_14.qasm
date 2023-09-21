OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(-0.13719288) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(4.5448137) q[1];
sx q[1];
rz(9.9546976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7512902) q[0];
sx q[0];
rz(-1.4634702) q[0];
sx q[0];
rz(-1.1885378) q[0];
rz(1.9723188) q[2];
sx q[2];
rz(-1.2520773) q[2];
sx q[2];
rz(2.4674494) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4094028) q[1];
sx q[1];
rz(-0.70530546) q[1];
sx q[1];
rz(0.7440872) q[1];
rz(1.3210117) q[3];
sx q[3];
rz(-2.3134109) q[3];
sx q[3];
rz(-3.0299377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4522176) q[2];
sx q[2];
rz(-1.8415425) q[2];
sx q[2];
rz(-0.33660647) q[2];
rz(-1.6254788) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(-1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34823725) q[0];
sx q[0];
rz(-1.1084778) q[0];
sx q[0];
rz(3.120378) q[0];
rz(-1.9477828) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(2.3056727) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056136925) q[0];
sx q[0];
rz(-1.5888927) q[0];
sx q[0];
rz(1.7096814) q[0];
rz(-1.0216653) q[2];
sx q[2];
rz(-1.9252535) q[2];
sx q[2];
rz(-2.3159388) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1509717) q[1];
sx q[1];
rz(-0.70578209) q[1];
sx q[1];
rz(-1.9981995) q[1];
x q[2];
rz(2.4017176) q[3];
sx q[3];
rz(-1.2736819) q[3];
sx q[3];
rz(-1.5418996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8643643) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(-1.7956087) q[2];
rz(-0.35955444) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42831746) q[0];
sx q[0];
rz(-1.0773766) q[0];
sx q[0];
rz(1.0536449) q[0];
rz(1.2288278) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(-2.704481) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1431883) q[0];
sx q[0];
rz(-0.30448738) q[0];
sx q[0];
rz(1.2782774) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24361165) q[2];
sx q[2];
rz(-1.1738452) q[2];
sx q[2];
rz(1.5543907) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25698173) q[1];
sx q[1];
rz(-1.4323438) q[1];
sx q[1];
rz(-2.715766) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0655572) q[3];
sx q[3];
rz(-1.5177739) q[3];
sx q[3];
rz(3.0474636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1221216) q[2];
sx q[2];
rz(-2.3601668) q[2];
sx q[2];
rz(1.0220698) q[2];
rz(1.2381037) q[3];
sx q[3];
rz(-0.3823897) q[3];
sx q[3];
rz(-0.4195655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5220752) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(2.1602901) q[0];
rz(-3.006382) q[1];
sx q[1];
rz(-1.0842666) q[1];
sx q[1];
rz(-0.19128004) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080973074) q[0];
sx q[0];
rz(-1.4220211) q[0];
sx q[0];
rz(-0.087555126) q[0];
x q[1];
rz(2.3046266) q[2];
sx q[2];
rz(-1.776473) q[2];
sx q[2];
rz(2.2107746) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0335238) q[1];
sx q[1];
rz(-0.78447693) q[1];
sx q[1];
rz(-0.14524059) q[1];
x q[2];
rz(2.5883834) q[3];
sx q[3];
rz(-0.60713967) q[3];
sx q[3];
rz(-2.7943484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4613351) q[2];
sx q[2];
rz(-0.98538435) q[2];
sx q[2];
rz(1.0106687) q[2];
rz(0.7615532) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33870944) q[0];
sx q[0];
rz(-2.8864679) q[0];
sx q[0];
rz(2.5849735) q[0];
rz(-3.026475) q[1];
sx q[1];
rz(-1.3373673) q[1];
sx q[1];
rz(-2.1690878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2199729) q[0];
sx q[0];
rz(-0.12244206) q[0];
sx q[0];
rz(-0.58453154) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33072492) q[2];
sx q[2];
rz(-0.84825584) q[2];
sx q[2];
rz(-1.5028138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79341187) q[1];
sx q[1];
rz(-0.23985292) q[1];
sx q[1];
rz(-2.1972149) q[1];
x q[2];
rz(0.95485165) q[3];
sx q[3];
rz(-1.6007489) q[3];
sx q[3];
rz(2.0314465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3187023) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(1.548432) q[2];
rz(-1.7758153) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(2.1877066) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71972972) q[0];
sx q[0];
rz(-1.2635764) q[0];
sx q[0];
rz(1.7156037) q[0];
rz(-1.0643719) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(-0.37429601) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6573338) q[0];
sx q[0];
rz(-1.9801635) q[0];
sx q[0];
rz(1.2249468) q[0];
rz(0.49403814) q[2];
sx q[2];
rz(-1.8015773) q[2];
sx q[2];
rz(-1.5766694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.53257886) q[1];
sx q[1];
rz(-0.62090579) q[1];
sx q[1];
rz(-2.6780374) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4738594) q[3];
sx q[3];
rz(-0.99620921) q[3];
sx q[3];
rz(2.7001065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2465683) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(-0.95823112) q[2];
rz(2.9124177) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(-0.62098256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.7763057) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(-0.90674415) q[0];
rz(-1.0892185) q[1];
sx q[1];
rz(-1.6420495) q[1];
sx q[1];
rz(-1.3100756) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.487264) q[0];
sx q[0];
rz(-1.1833106) q[0];
sx q[0];
rz(1.6261473) q[0];
rz(1.9016978) q[2];
sx q[2];
rz(-2.2563997) q[2];
sx q[2];
rz(-0.14724018) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20565614) q[1];
sx q[1];
rz(-1.9034791) q[1];
sx q[1];
rz(-1.7722305) q[1];
rz(-1.5433611) q[3];
sx q[3];
rz(-1.0155639) q[3];
sx q[3];
rz(-0.96019402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92581302) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(1.4833935) q[2];
rz(0.27967927) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(-3.083995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9782372) q[0];
sx q[0];
rz(-1.5196479) q[0];
sx q[0];
rz(-0.21959198) q[0];
rz(-2.638468) q[1];
sx q[1];
rz(-0.88880912) q[1];
sx q[1];
rz(2.2917152) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06892878) q[0];
sx q[0];
rz(-1.7964296) q[0];
sx q[0];
rz(0.16107852) q[0];
rz(-0.093896534) q[2];
sx q[2];
rz(-2.2517423) q[2];
sx q[2];
rz(2.4927793) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5819468) q[1];
sx q[1];
rz(-1.6053182) q[1];
sx q[1];
rz(-0.80969651) q[1];
rz(-pi) q[2];
rz(-3.0181846) q[3];
sx q[3];
rz(-1.3128237) q[3];
sx q[3];
rz(-2.3464399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5510817) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(1.760651) q[2];
rz(-0.75602174) q[3];
sx q[3];
rz(-0.20320007) q[3];
sx q[3];
rz(-0.35593629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6324156) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(-0.40503043) q[0];
rz(-2.6889154) q[1];
sx q[1];
rz(-0.98222268) q[1];
sx q[1];
rz(1.8639494) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4687846) q[0];
sx q[0];
rz(-2.0622258) q[0];
sx q[0];
rz(1.8322893) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69127609) q[2];
sx q[2];
rz(-1.1862159) q[2];
sx q[2];
rz(2.861475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.35754044) q[1];
sx q[1];
rz(-1.5338384) q[1];
sx q[1];
rz(2.0164911) q[1];
rz(2.7550335) q[3];
sx q[3];
rz(-2.2694526) q[3];
sx q[3];
rz(-2.3545635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8032802) q[2];
sx q[2];
rz(-0.40643224) q[2];
sx q[2];
rz(2.7837616) q[2];
rz(1.7221649) q[3];
sx q[3];
rz(-1.8678886) q[3];
sx q[3];
rz(1.0740124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91838592) q[0];
sx q[0];
rz(-3.0637488) q[0];
sx q[0];
rz(0.11225587) q[0];
rz(2.2414801) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(0.21044883) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8716547) q[0];
sx q[0];
rz(-2.6007814) q[0];
sx q[0];
rz(2.7842115) q[0];
rz(-1.2082645) q[2];
sx q[2];
rz(-0.4368442) q[2];
sx q[2];
rz(-0.7561965) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.12293832) q[1];
sx q[1];
rz(-1.8784874) q[1];
sx q[1];
rz(-1.3535181) q[1];
rz(-pi) q[2];
rz(-2.8912192) q[3];
sx q[3];
rz(-1.3406521) q[3];
sx q[3];
rz(-1.9104513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.18008733) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(-1.5853184) q[2];
rz(-1.2735584) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(2.9343228) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5719941) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(-2.3251484) q[1];
sx q[1];
rz(-1.888231) q[1];
sx q[1];
rz(2.9838557) q[1];
rz(2.0408761) q[2];
sx q[2];
rz(-1.4224978) q[2];
sx q[2];
rz(-2.9329621) q[2];
rz(0.3420842) q[3];
sx q[3];
rz(-0.87089201) q[3];
sx q[3];
rz(1.0637829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];