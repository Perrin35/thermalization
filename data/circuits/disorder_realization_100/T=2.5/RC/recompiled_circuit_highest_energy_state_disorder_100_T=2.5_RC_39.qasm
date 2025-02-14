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
rz(-0.21021357) q[0];
sx q[0];
rz(-2.7855594) q[0];
sx q[0];
rz(1.5176679) q[0];
rz(-1.3104982) q[1];
sx q[1];
rz(-0.85135353) q[1];
sx q[1];
rz(0.89827615) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0097281) q[0];
sx q[0];
rz(-0.72355958) q[0];
sx q[0];
rz(-1.0148763) q[0];
rz(0.29062985) q[2];
sx q[2];
rz(-1.0950452) q[2];
sx q[2];
rz(-1.0889183) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9834302) q[1];
sx q[1];
rz(-1.753141) q[1];
sx q[1];
rz(-1.4901864) q[1];
rz(-pi) q[2];
rz(-0.0086011767) q[3];
sx q[3];
rz(-1.4066753) q[3];
sx q[3];
rz(-0.65879956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6015168) q[2];
sx q[2];
rz(-2.8105152) q[2];
sx q[2];
rz(-0.69851056) q[2];
rz(1.6554333) q[3];
sx q[3];
rz(-0.58659068) q[3];
sx q[3];
rz(1.9894039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.0450714) q[0];
sx q[0];
rz(-0.69913816) q[0];
sx q[0];
rz(-2.844098) q[0];
rz(2.7104764) q[1];
sx q[1];
rz(-2.2747206) q[1];
sx q[1];
rz(-1.3131622) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8968412) q[0];
sx q[0];
rz(-0.8936106) q[0];
sx q[0];
rz(2.0170101) q[0];
rz(-1.976491) q[2];
sx q[2];
rz(-1.5303474) q[2];
sx q[2];
rz(1.2834398) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0894236) q[1];
sx q[1];
rz(-2.3681297) q[1];
sx q[1];
rz(-1.5599584) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24784577) q[3];
sx q[3];
rz(-1.447926) q[3];
sx q[3];
rz(2.7932699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7737274) q[2];
sx q[2];
rz(-0.44168681) q[2];
sx q[2];
rz(2.3750677) q[2];
rz(0.68486989) q[3];
sx q[3];
rz(-1.6133512) q[3];
sx q[3];
rz(-0.49055704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8399452) q[0];
sx q[0];
rz(-1.8122346) q[0];
sx q[0];
rz(-2.8369821) q[0];
rz(-2.1412762) q[1];
sx q[1];
rz(-1.1023003) q[1];
sx q[1];
rz(2.2432378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2294126) q[0];
sx q[0];
rz(-2.1026975) q[0];
sx q[0];
rz(0.096166111) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86849536) q[2];
sx q[2];
rz(-1.8284594) q[2];
sx q[2];
rz(2.9106377) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5872941) q[1];
sx q[1];
rz(-2.0600494) q[1];
sx q[1];
rz(2.4913408) q[1];
rz(-pi) q[2];
rz(1.0482084) q[3];
sx q[3];
rz(-1.6500123) q[3];
sx q[3];
rz(-1.0635215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0072713) q[2];
sx q[2];
rz(-1.4192702) q[2];
sx q[2];
rz(-2.814494) q[2];
rz(-1.6359811) q[3];
sx q[3];
rz(-1.4666731) q[3];
sx q[3];
rz(-1.8065642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.026767749) q[0];
sx q[0];
rz(-0.46634316) q[0];
sx q[0];
rz(-0.53792167) q[0];
rz(2.3665358) q[1];
sx q[1];
rz(-1.8102976) q[1];
sx q[1];
rz(0.34781003) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3191158) q[0];
sx q[0];
rz(-2.4317784) q[0];
sx q[0];
rz(0.7546203) q[0];
rz(3.0362282) q[2];
sx q[2];
rz(-1.67571) q[2];
sx q[2];
rz(2.9199469) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0671652) q[1];
sx q[1];
rz(-0.8377155) q[1];
sx q[1];
rz(-0.58372402) q[1];
x q[2];
rz(-1.9815481) q[3];
sx q[3];
rz(-2.7166945) q[3];
sx q[3];
rz(-1.9187669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20370087) q[2];
sx q[2];
rz(-0.72526473) q[2];
sx q[2];
rz(2.4791278) q[2];
rz(2.0857816) q[3];
sx q[3];
rz(-0.31848389) q[3];
sx q[3];
rz(-1.2663579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98589677) q[0];
sx q[0];
rz(-0.53845423) q[0];
sx q[0];
rz(1.023531) q[0];
rz(-1.0515949) q[1];
sx q[1];
rz(-1.4902427) q[1];
sx q[1];
rz(-0.016062707) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9413105) q[0];
sx q[0];
rz(-0.56277983) q[0];
sx q[0];
rz(-1.3154696) q[0];
rz(0.93463411) q[2];
sx q[2];
rz(-2.155288) q[2];
sx q[2];
rz(-0.57022695) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6391568) q[1];
sx q[1];
rz(-1.0519439) q[1];
sx q[1];
rz(-1.3154037) q[1];
x q[2];
rz(3.06918) q[3];
sx q[3];
rz(-1.2999897) q[3];
sx q[3];
rz(-0.31596247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.84231275) q[2];
sx q[2];
rz(-0.88025847) q[2];
sx q[2];
rz(-0.24097815) q[2];
rz(1.109451) q[3];
sx q[3];
rz(-1.2883319) q[3];
sx q[3];
rz(-0.39458767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1269162) q[0];
sx q[0];
rz(-0.83634818) q[0];
sx q[0];
rz(-0.29689223) q[0];
rz(-1.1709921) q[1];
sx q[1];
rz(-1.820727) q[1];
sx q[1];
rz(2.6062633) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2123665) q[0];
sx q[0];
rz(-1.3844212) q[0];
sx q[0];
rz(1.7586238) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8814323) q[2];
sx q[2];
rz(-2.363003) q[2];
sx q[2];
rz(-1.006135) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3228906) q[1];
sx q[1];
rz(-1.7826347) q[1];
sx q[1];
rz(1.0781652) q[1];
rz(-pi) q[2];
rz(3.0020079) q[3];
sx q[3];
rz(-1.4142766) q[3];
sx q[3];
rz(-0.28790441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8239596) q[2];
sx q[2];
rz(-2.3667658) q[2];
sx q[2];
rz(-0.63014692) q[2];
rz(-1.0050425) q[3];
sx q[3];
rz(-2.9229087) q[3];
sx q[3];
rz(2.4439243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.180535) q[0];
sx q[0];
rz(-3.0735885) q[0];
sx q[0];
rz(1.8649752) q[0];
rz(-2.0969773) q[1];
sx q[1];
rz(-1.4199384) q[1];
sx q[1];
rz(2.9081664) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7395494) q[0];
sx q[0];
rz(-1.6130216) q[0];
sx q[0];
rz(2.2916138) q[0];
rz(2.7415375) q[2];
sx q[2];
rz(-2.4597485) q[2];
sx q[2];
rz(1.2712196) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0580768) q[1];
sx q[1];
rz(-2.3928071) q[1];
sx q[1];
rz(-0.056549055) q[1];
rz(-3.0393321) q[3];
sx q[3];
rz(-0.58057154) q[3];
sx q[3];
rz(1.5619265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.476568) q[2];
sx q[2];
rz(-1.2580405) q[2];
sx q[2];
rz(-0.27102077) q[2];
rz(-2.792231) q[3];
sx q[3];
rz(-0.91785279) q[3];
sx q[3];
rz(-2.5640986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2006328) q[0];
sx q[0];
rz(-2.2043493) q[0];
sx q[0];
rz(-1.3295133) q[0];
rz(-0.26893523) q[1];
sx q[1];
rz(-0.66489282) q[1];
sx q[1];
rz(1.7609133) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5047245) q[0];
sx q[0];
rz(-0.32384017) q[0];
sx q[0];
rz(2.25405) q[0];
rz(-pi) q[1];
rz(3.1166638) q[2];
sx q[2];
rz(-2.5453794) q[2];
sx q[2];
rz(0.91170694) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.096752874) q[1];
sx q[1];
rz(-2.1659454) q[1];
sx q[1];
rz(0.19697856) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56048067) q[3];
sx q[3];
rz(-1.4541601) q[3];
sx q[3];
rz(2.9612535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55714503) q[2];
sx q[2];
rz(-2.2717805) q[2];
sx q[2];
rz(2.4073041) q[2];
rz(2.0860489) q[3];
sx q[3];
rz(-1.5493834) q[3];
sx q[3];
rz(-0.9744823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.664481) q[0];
sx q[0];
rz(-2.8900914) q[0];
sx q[0];
rz(2.2736881) q[0];
rz(-1.3033298) q[1];
sx q[1];
rz(-1.5497327) q[1];
sx q[1];
rz(-0.23095362) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8036007) q[0];
sx q[0];
rz(-0.70728318) q[0];
sx q[0];
rz(2.4770205) q[0];
rz(-2.7074293) q[2];
sx q[2];
rz(-1.4462398) q[2];
sx q[2];
rz(-1.9366154) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8924397) q[1];
sx q[1];
rz(-2.0241535) q[1];
sx q[1];
rz(-1.3282177) q[1];
x q[2];
rz(-0.001570985) q[3];
sx q[3];
rz(-0.6720619) q[3];
sx q[3];
rz(1.5379958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2086198) q[2];
sx q[2];
rz(-1.695637) q[2];
sx q[2];
rz(0.43819532) q[2];
rz(2.4402601) q[3];
sx q[3];
rz(-2.0739136) q[3];
sx q[3];
rz(-0.35593885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0375131) q[0];
sx q[0];
rz(-0.47261819) q[0];
sx q[0];
rz(-1.0802826) q[0];
rz(-1.0516306) q[1];
sx q[1];
rz(-2.0104505) q[1];
sx q[1];
rz(-2.0463321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29275492) q[0];
sx q[0];
rz(-1.0612604) q[0];
sx q[0];
rz(-1.900702) q[0];
rz(-0.3682179) q[2];
sx q[2];
rz(-1.5309586) q[2];
sx q[2];
rz(1.4152272) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0174696) q[1];
sx q[1];
rz(-0.7354722) q[1];
sx q[1];
rz(-1.7361438) q[1];
rz(2.0729823) q[3];
sx q[3];
rz(-2.1081703) q[3];
sx q[3];
rz(0.93056941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.1209391) q[2];
sx q[2];
rz(-1.8755269) q[2];
sx q[2];
rz(0.90547639) q[2];
rz(-0.71546537) q[3];
sx q[3];
rz(-2.4162636) q[3];
sx q[3];
rz(0.65413094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.41659551) q[0];
sx q[0];
rz(-1.7323957) q[0];
sx q[0];
rz(0.61687627) q[0];
rz(2.0116518) q[1];
sx q[1];
rz(-1.8604953) q[1];
sx q[1];
rz(-1.4166191) q[1];
rz(-2.5196709) q[2];
sx q[2];
rz(-2.1991232) q[2];
sx q[2];
rz(-0.029673619) q[2];
rz(-2.1917584) q[3];
sx q[3];
rz(-0.25706188) q[3];
sx q[3];
rz(1.2973447) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
