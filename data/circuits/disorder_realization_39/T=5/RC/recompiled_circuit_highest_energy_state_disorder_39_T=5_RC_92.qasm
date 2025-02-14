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
rz(0.14035913) q[0];
sx q[0];
rz(4.601534) q[0];
sx q[0];
rz(10.277716) q[0];
rz(-1.1879022) q[1];
sx q[1];
rz(-0.24628425) q[1];
sx q[1];
rz(-2.8077717) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.649619) q[0];
sx q[0];
rz(-2.5707757) q[0];
sx q[0];
rz(1.094857) q[0];
x q[1];
rz(-0.63392459) q[2];
sx q[2];
rz(-1.9854387) q[2];
sx q[2];
rz(-0.3699322) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8454518) q[1];
sx q[1];
rz(-2.9036464) q[1];
sx q[1];
rz(0.68745698) q[1];
rz(-pi) q[2];
rz(-3.0404011) q[3];
sx q[3];
rz(-1.4863761) q[3];
sx q[3];
rz(-2.3197378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21981257) q[2];
sx q[2];
rz(-2.2647936) q[2];
sx q[2];
rz(0.63966695) q[2];
rz(-0.93846792) q[3];
sx q[3];
rz(-0.42249051) q[3];
sx q[3];
rz(1.2708906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.4939782) q[0];
sx q[0];
rz(-3.0632126) q[0];
sx q[0];
rz(1.6139503) q[0];
rz(-2.899462) q[1];
sx q[1];
rz(-2.1277728) q[1];
sx q[1];
rz(2.4193144) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9856208) q[0];
sx q[0];
rz(-1.9160144) q[0];
sx q[0];
rz(1.9916608) q[0];
rz(-pi) q[1];
rz(0.26210384) q[2];
sx q[2];
rz(-2.3293709) q[2];
sx q[2];
rz(-1.1103528) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5879155) q[1];
sx q[1];
rz(-2.0688631) q[1];
sx q[1];
rz(0.95446511) q[1];
rz(2.9728209) q[3];
sx q[3];
rz(-0.95484551) q[3];
sx q[3];
rz(1.7273161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0904514) q[2];
sx q[2];
rz(-2.800056) q[2];
sx q[2];
rz(-0.086645834) q[2];
rz(-2.5610949) q[3];
sx q[3];
rz(-1.2248421) q[3];
sx q[3];
rz(-2.1897924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8101623) q[0];
sx q[0];
rz(-2.0482752) q[0];
sx q[0];
rz(1.4587559) q[0];
rz(3.0259865) q[1];
sx q[1];
rz(-1.9435725) q[1];
sx q[1];
rz(-2.9748532) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34761679) q[0];
sx q[0];
rz(-2.1361094) q[0];
sx q[0];
rz(0.24576776) q[0];
rz(-pi) q[1];
rz(0.42129321) q[2];
sx q[2];
rz(-1.2308106) q[2];
sx q[2];
rz(0.30145633) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3003242) q[1];
sx q[1];
rz(-2.1144697) q[1];
sx q[1];
rz(-0.44692301) q[1];
rz(-1.429168) q[3];
sx q[3];
rz(-2.3109791) q[3];
sx q[3];
rz(0.52799559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.705767) q[2];
sx q[2];
rz(-0.21958084) q[2];
sx q[2];
rz(-1.1337918) q[2];
rz(0.052915834) q[3];
sx q[3];
rz(-1.2727126) q[3];
sx q[3];
rz(-2.4165966) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13919203) q[0];
sx q[0];
rz(-0.63183689) q[0];
sx q[0];
rz(-2.8644526) q[0];
rz(0.78284043) q[1];
sx q[1];
rz(-1.6354086) q[1];
sx q[1];
rz(-2.7908819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6576516) q[0];
sx q[0];
rz(-0.707905) q[0];
sx q[0];
rz(-1.4483676) q[0];
rz(1.6007623) q[2];
sx q[2];
rz(-1.2453516) q[2];
sx q[2];
rz(-1.1496227) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8102032) q[1];
sx q[1];
rz(-2.1955829) q[1];
sx q[1];
rz(0.58831711) q[1];
rz(-pi) q[2];
rz(1.1923741) q[3];
sx q[3];
rz(-3.0251092) q[3];
sx q[3];
rz(2.0476215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1002525) q[2];
sx q[2];
rz(-1.2790054) q[2];
sx q[2];
rz(-1.1502385) q[2];
rz(-2.9669115) q[3];
sx q[3];
rz(-2.8476069) q[3];
sx q[3];
rz(3.0512419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-3.1133076) q[0];
sx q[0];
rz(-2.568013) q[0];
sx q[0];
rz(1.9453402) q[0];
rz(-0.47736564) q[1];
sx q[1];
rz(-2.3148675) q[1];
sx q[1];
rz(-0.54944077) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4081037) q[0];
sx q[0];
rz(-2.2814352) q[0];
sx q[0];
rz(-2.3657794) q[0];
rz(-pi) q[1];
rz(-2.7844593) q[2];
sx q[2];
rz(-2.5741842) q[2];
sx q[2];
rz(1.0536989) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3122299) q[1];
sx q[1];
rz(-0.42604241) q[1];
sx q[1];
rz(1.8660353) q[1];
x q[2];
rz(-2.0277689) q[3];
sx q[3];
rz(-1.1871871) q[3];
sx q[3];
rz(-0.76573223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4598733) q[2];
sx q[2];
rz(-0.37178603) q[2];
sx q[2];
rz(0.2447153) q[2];
rz(-1.6230445) q[3];
sx q[3];
rz(-1.1145498) q[3];
sx q[3];
rz(-3.0486619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2731584) q[0];
sx q[0];
rz(-2.1491829) q[0];
sx q[0];
rz(-0.42700818) q[0];
rz(1.931841) q[1];
sx q[1];
rz(-1.6170343) q[1];
sx q[1];
rz(0.0078113656) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2586827) q[0];
sx q[0];
rz(-2.3091303) q[0];
sx q[0];
rz(-0.29405221) q[0];
rz(-pi) q[1];
rz(1.5862238) q[2];
sx q[2];
rz(-2.5944105) q[2];
sx q[2];
rz(1.6268961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5558034) q[1];
sx q[1];
rz(-0.92098707) q[1];
sx q[1];
rz(1.5073677) q[1];
x q[2];
rz(2.8704371) q[3];
sx q[3];
rz(-0.89224766) q[3];
sx q[3];
rz(-0.0038879768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.740364) q[2];
sx q[2];
rz(-0.87330356) q[2];
sx q[2];
rz(-2.002423) q[2];
rz(0.27979699) q[3];
sx q[3];
rz(-1.8179025) q[3];
sx q[3];
rz(-1.4626224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6497659) q[0];
sx q[0];
rz(-0.38182807) q[0];
sx q[0];
rz(2.5873798) q[0];
rz(0.062049374) q[1];
sx q[1];
rz(-0.65826145) q[1];
sx q[1];
rz(2.2668692) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.301247) q[0];
sx q[0];
rz(-0.78784993) q[0];
sx q[0];
rz(-1.3993174) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.563233) q[2];
sx q[2];
rz(-1.5652135) q[2];
sx q[2];
rz(-1.1712318) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3530423) q[1];
sx q[1];
rz(-1.0354831) q[1];
sx q[1];
rz(-2.8393406) q[1];
rz(-pi) q[2];
rz(-2.234972) q[3];
sx q[3];
rz(-1.2620838) q[3];
sx q[3];
rz(1.3332092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41500652) q[2];
sx q[2];
rz(-1.0174624) q[2];
sx q[2];
rz(0.93878186) q[2];
rz(-0.69432652) q[3];
sx q[3];
rz(-0.66803473) q[3];
sx q[3];
rz(2.9850849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.071455) q[0];
sx q[0];
rz(-2.5818765) q[0];
sx q[0];
rz(-2.8330084) q[0];
rz(1.6584819) q[1];
sx q[1];
rz(-1.7959271) q[1];
sx q[1];
rz(-2.9187091) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4083129) q[0];
sx q[0];
rz(-2.2594497) q[0];
sx q[0];
rz(1.0802425) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7168458) q[2];
sx q[2];
rz(-0.24492376) q[2];
sx q[2];
rz(0.67451678) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0012344) q[1];
sx q[1];
rz(-1.8158633) q[1];
sx q[1];
rz(3.1015741) q[1];
rz(-pi) q[2];
rz(2.4588575) q[3];
sx q[3];
rz(-2.020524) q[3];
sx q[3];
rz(2.3537824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3234723) q[2];
sx q[2];
rz(-1.0264531) q[2];
sx q[2];
rz(0.66413122) q[2];
rz(-1.9984261) q[3];
sx q[3];
rz(-1.6616471) q[3];
sx q[3];
rz(-2.0460879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42300647) q[0];
sx q[0];
rz(-1.9976595) q[0];
sx q[0];
rz(-0.018420694) q[0];
rz(2.9711235) q[1];
sx q[1];
rz(-1.2455995) q[1];
sx q[1];
rz(-2.9439994) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22329837) q[0];
sx q[0];
rz(-2.9726148) q[0];
sx q[0];
rz(-1.5150945) q[0];
rz(-pi) q[1];
rz(-1.204448) q[2];
sx q[2];
rz(-1.5218668) q[2];
sx q[2];
rz(1.4885224) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9925633) q[1];
sx q[1];
rz(-0.32900336) q[1];
sx q[1];
rz(2.002153) q[1];
x q[2];
rz(-0.13473265) q[3];
sx q[3];
rz(-2.7425457) q[3];
sx q[3];
rz(1.4302554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1175055) q[2];
sx q[2];
rz(-1.0988289) q[2];
sx q[2];
rz(-2.7433024) q[2];
rz(2.6309218) q[3];
sx q[3];
rz(-1.4887678) q[3];
sx q[3];
rz(-2.2834856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2354105) q[0];
sx q[0];
rz(-2.372083) q[0];
sx q[0];
rz(0.19943516) q[0];
rz(-0.27345744) q[1];
sx q[1];
rz(-2.7077935) q[1];
sx q[1];
rz(-1.010703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.511925) q[0];
sx q[0];
rz(-2.2154847) q[0];
sx q[0];
rz(1.900338) q[0];
x q[1];
rz(1.550714) q[2];
sx q[2];
rz(-1.3354567) q[2];
sx q[2];
rz(-2.2724336) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6335754) q[1];
sx q[1];
rz(-2.4035638) q[1];
sx q[1];
rz(-0.077390093) q[1];
rz(1.8685249) q[3];
sx q[3];
rz(-1.862251) q[3];
sx q[3];
rz(3.1133661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.79591862) q[2];
sx q[2];
rz(-0.60563874) q[2];
sx q[2];
rz(-2.3460713) q[2];
rz(-2.7703721) q[3];
sx q[3];
rz(-1.755654) q[3];
sx q[3];
rz(2.8871239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026074792) q[0];
sx q[0];
rz(-1.4160897) q[0];
sx q[0];
rz(-1.8427451) q[0];
rz(-0.55116354) q[1];
sx q[1];
rz(-1.9438585) q[1];
sx q[1];
rz(1.9895947) q[1];
rz(0.25210512) q[2];
sx q[2];
rz(-0.74130836) q[2];
sx q[2];
rz(3.1234968) q[2];
rz(-2.9007111) q[3];
sx q[3];
rz(-1.1987562) q[3];
sx q[3];
rz(-0.22507122) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
