OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8496348) q[0];
sx q[0];
rz(-0.44009527) q[0];
sx q[0];
rz(0.13719288) q[0];
rz(-1.7358915) q[1];
sx q[1];
rz(-1.403221) q[1];
sx q[1];
rz(-0.52991968) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0041381) q[0];
sx q[0];
rz(-1.9507427) q[0];
sx q[0];
rz(-3.0259892) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.797384) q[2];
sx q[2];
rz(-1.1905626) q[2];
sx q[2];
rz(2.3772079) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4094028) q[1];
sx q[1];
rz(-2.4362872) q[1];
sx q[1];
rz(-2.3975055) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3210117) q[3];
sx q[3];
rz(-0.82818177) q[3];
sx q[3];
rz(3.0299377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4522176) q[2];
sx q[2];
rz(-1.8415425) q[2];
sx q[2];
rz(-2.8049862) q[2];
rz(-1.5161139) q[3];
sx q[3];
rz(-2.5879526) q[3];
sx q[3];
rz(1.6158993) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7933554) q[0];
sx q[0];
rz(-1.1084778) q[0];
sx q[0];
rz(-0.021214699) q[0];
rz(-1.1938098) q[1];
sx q[1];
rz(-1.0394916) q[1];
sx q[1];
rz(-0.83591998) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3859235) q[0];
sx q[0];
rz(-3.0015411) q[0];
sx q[0];
rz(-1.4408017) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0216653) q[2];
sx q[2];
rz(-1.9252535) q[2];
sx q[2];
rz(2.3159388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6902496) q[1];
sx q[1];
rz(-2.2022044) q[1];
sx q[1];
rz(-0.33957014) q[1];
rz(-0.73987506) q[3];
sx q[3];
rz(-1.8679108) q[3];
sx q[3];
rz(-1.5996931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2772284) q[2];
sx q[2];
rz(-1.1479062) q[2];
sx q[2];
rz(1.345984) q[2];
rz(0.35955444) q[3];
sx q[3];
rz(-0.94272009) q[3];
sx q[3];
rz(0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7132752) q[0];
sx q[0];
rz(-1.0773766) q[0];
sx q[0];
rz(2.0879478) q[0];
rz(1.9127649) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(-0.4371117) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99840435) q[0];
sx q[0];
rz(-2.8371053) q[0];
sx q[0];
rz(-1.2782774) q[0];
rz(-1.0486629) q[2];
sx q[2];
rz(-2.6792567) q[2];
sx q[2];
rz(-2.1585652) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8846109) q[1];
sx q[1];
rz(-1.4323438) q[1];
sx q[1];
rz(-2.715766) q[1];
rz(1.6821074) q[3];
sx q[3];
rz(-2.6442332) q[3];
sx q[3];
rz(1.5670083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.019471021) q[2];
sx q[2];
rz(-2.3601668) q[2];
sx q[2];
rz(1.0220698) q[2];
rz(1.9034889) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(-0.4195655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6195174) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(-0.98130256) q[0];
rz(-0.13521067) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(-0.19128004) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45368886) q[0];
sx q[0];
rz(-0.17246248) q[0];
sx q[0];
rz(2.0989291) q[0];
rz(2.3046266) q[2];
sx q[2];
rz(-1.3651197) q[2];
sx q[2];
rz(-2.2107746) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7820218) q[1];
sx q[1];
rz(-1.6732209) q[1];
sx q[1];
rz(-2.3624079) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55320923) q[3];
sx q[3];
rz(-0.60713967) q[3];
sx q[3];
rz(0.34724423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4613351) q[2];
sx q[2];
rz(-0.98538435) q[2];
sx q[2];
rz(-1.0106687) q[2];
rz(-2.3800395) q[3];
sx q[3];
rz(-1.9617617) q[3];
sx q[3];
rz(-0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8028832) q[0];
sx q[0];
rz(-2.8864679) q[0];
sx q[0];
rz(0.55661911) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(0.97250485) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33362493) q[0];
sx q[0];
rz(-1.6728314) q[0];
sx q[0];
rz(1.6385965) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3210578) q[2];
sx q[2];
rz(-1.3247326) q[2];
sx q[2];
rz(2.8503502) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.433686) q[1];
sx q[1];
rz(-1.3771332) q[1];
sx q[1];
rz(-0.14240264) q[1];
rz(1.6226107) q[3];
sx q[3];
rz(-0.6165781) q[3];
sx q[3];
rz(-2.6386564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.82289034) q[2];
sx q[2];
rz(-2.0501142) q[2];
sx q[2];
rz(1.5931607) q[2];
rz(-1.7758153) q[3];
sx q[3];
rz(-0.32309353) q[3];
sx q[3];
rz(0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4218629) q[0];
sx q[0];
rz(-1.2635764) q[0];
sx q[0];
rz(1.4259889) q[0];
rz(2.0772207) q[1];
sx q[1];
rz(-1.0168889) q[1];
sx q[1];
rz(-2.7672966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4842589) q[0];
sx q[0];
rz(-1.1614292) q[0];
sx q[0];
rz(1.9166458) q[0];
rz(-pi) q[1];
rz(1.8315973) q[2];
sx q[2];
rz(-2.0506095) q[2];
sx q[2];
rz(3.0248883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7172076) q[1];
sx q[1];
rz(-1.8339515) q[1];
sx q[1];
rz(2.5724263) q[1];
x q[2];
rz(-2.3338942) q[3];
sx q[3];
rz(-2.2904615) q[3];
sx q[3];
rz(-2.6157275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8950243) q[2];
sx q[2];
rz(-2.4763069) q[2];
sx q[2];
rz(-0.95823112) q[2];
rz(2.9124177) q[3];
sx q[3];
rz(-1.6848247) q[3];
sx q[3];
rz(-2.5206101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7763057) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(0.90674415) q[0];
rz(-2.0523741) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(-1.3100756) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6543286) q[0];
sx q[0];
rz(-1.9582821) q[0];
sx q[0];
rz(-1.5154454) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9016978) q[2];
sx q[2];
rz(-2.2563997) q[2];
sx q[2];
rz(-2.9943525) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9359365) q[1];
sx q[1];
rz(-1.9034791) q[1];
sx q[1];
rz(-1.7722305) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[1];
rz(-0.92581302) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(1.4833935) q[2];
rz(-2.8619134) q[3];
sx q[3];
rz(-2.1488583) q[3];
sx q[3];
rz(3.083995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9782372) q[0];
sx q[0];
rz(-1.6219448) q[0];
sx q[0];
rz(0.21959198) q[0];
rz(-2.638468) q[1];
sx q[1];
rz(-2.2527835) q[1];
sx q[1];
rz(-2.2917152) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06892878) q[0];
sx q[0];
rz(-1.7964296) q[0];
sx q[0];
rz(0.16107852) q[0];
x q[1];
rz(-2.2539027) q[2];
sx q[2];
rz(-1.6437093) q[2];
sx q[2];
rz(2.2788252) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0440158) q[1];
sx q[1];
rz(-2.3313287) q[1];
sx q[1];
rz(-3.0939328) q[1];
rz(-pi) q[2];
rz(0.12340801) q[3];
sx q[3];
rz(-1.3128237) q[3];
sx q[3];
rz(-2.3464399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59051096) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(1.3809416) q[2];
rz(-0.75602174) q[3];
sx q[3];
rz(-0.20320007) q[3];
sx q[3];
rz(-0.35593629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6324156) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(0.40503043) q[0];
rz(0.45267725) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(1.2776432) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4687846) q[0];
sx q[0];
rz(-1.0793669) q[0];
sx q[0];
rz(1.3093033) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69127609) q[2];
sx q[2];
rz(-1.1862159) q[2];
sx q[2];
rz(0.28011766) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9459878) q[1];
sx q[1];
rz(-1.1254278) q[1];
sx q[1];
rz(-0.040954879) q[1];
rz(-2.7550335) q[3];
sx q[3];
rz(-2.2694526) q[3];
sx q[3];
rz(2.3545635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3383125) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(0.35783106) q[2];
rz(-1.7221649) q[3];
sx q[3];
rz(-1.8678886) q[3];
sx q[3];
rz(2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2232067) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(-3.0293368) q[0];
rz(2.2414801) q[1];
sx q[1];
rz(-2.0745514) q[1];
sx q[1];
rz(2.9311438) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99104098) q[0];
sx q[0];
rz(-1.3897087) q[0];
sx q[0];
rz(-0.51245706) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9333282) q[2];
sx q[2];
rz(-2.7047485) q[2];
sx q[2];
rz(-0.7561965) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0186543) q[1];
sx q[1];
rz(-1.2631053) q[1];
sx q[1];
rz(1.3535181) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3335161) q[3];
sx q[3];
rz(-1.3271601) q[3];
sx q[3];
rz(0.39792774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9615053) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(1.5562742) q[2];
rz(-1.8680343) q[3];
sx q[3];
rz(-2.5189416) q[3];
sx q[3];
rz(-0.20726985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5719941) q[0];
sx q[0];
rz(-2.270569) q[0];
sx q[0];
rz(1.7763174) q[0];
rz(0.81644425) q[1];
sx q[1];
rz(-1.888231) q[1];
sx q[1];
rz(2.9838557) q[1];
rz(-2.9755637) q[2];
sx q[2];
rz(-2.0353073) q[2];
sx q[2];
rz(-1.4370949) q[2];
rz(2.7995085) q[3];
sx q[3];
rz(-2.2707006) q[3];
sx q[3];
rz(-2.0778098) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
