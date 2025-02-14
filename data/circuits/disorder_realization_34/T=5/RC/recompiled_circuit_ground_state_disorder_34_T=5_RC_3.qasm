OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7393957) q[0];
sx q[0];
rz(4.0255044) q[0];
sx q[0];
rz(8.5691353) q[0];
rz(-0.17172509) q[1];
sx q[1];
rz(-0.11556927) q[1];
sx q[1];
rz(2.5791383) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64096686) q[0];
sx q[0];
rz(-1.62407) q[0];
sx q[0];
rz(1.9356273) q[0];
x q[1];
rz(-2.0777656) q[2];
sx q[2];
rz(-1.7071144) q[2];
sx q[2];
rz(0.16198128) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2515638) q[1];
sx q[1];
rz(-1.152413) q[1];
sx q[1];
rz(2.0189925) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1645557) q[3];
sx q[3];
rz(-1.2066168) q[3];
sx q[3];
rz(1.6842357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2354551) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(2.3423024) q[2];
rz(2.6702787) q[3];
sx q[3];
rz(-2.2282232) q[3];
sx q[3];
rz(-2.2057064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1021295) q[0];
sx q[0];
rz(-1.3251745) q[0];
sx q[0];
rz(1.7934196) q[0];
rz(-1.3735636) q[1];
sx q[1];
rz(-1.1496239) q[1];
sx q[1];
rz(-1.9893533) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3593529) q[0];
sx q[0];
rz(-0.77762654) q[0];
sx q[0];
rz(-1.9944784) q[0];
rz(-pi) q[1];
rz(2.4068314) q[2];
sx q[2];
rz(-0.82345357) q[2];
sx q[2];
rz(-2.6121989) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5630398) q[1];
sx q[1];
rz(-0.7078979) q[1];
sx q[1];
rz(2.8824174) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1227319) q[3];
sx q[3];
rz(-2.0235007) q[3];
sx q[3];
rz(2.798003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.79933244) q[2];
sx q[2];
rz(-2.7228184) q[2];
sx q[2];
rz(2.2104134) q[2];
rz(-0.10654199) q[3];
sx q[3];
rz(-1.1139161) q[3];
sx q[3];
rz(2.54134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057864144) q[0];
sx q[0];
rz(-1.7513542) q[0];
sx q[0];
rz(2.6237543) q[0];
rz(2.2528516) q[1];
sx q[1];
rz(-0.70777142) q[1];
sx q[1];
rz(-0.4471561) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6155375) q[0];
sx q[0];
rz(-1.3096022) q[0];
sx q[0];
rz(-1.636807) q[0];
x q[1];
rz(-2.4141623) q[2];
sx q[2];
rz(-1.1377678) q[2];
sx q[2];
rz(-0.43011452) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0886317) q[1];
sx q[1];
rz(-1.4354032) q[1];
sx q[1];
rz(-0.9779344) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6154061) q[3];
sx q[3];
rz(-1.8136214) q[3];
sx q[3];
rz(-0.015004166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37457028) q[2];
sx q[2];
rz(-2.3775103) q[2];
sx q[2];
rz(-2.6089597) q[2];
rz(0.24752188) q[3];
sx q[3];
rz(-2.4041924) q[3];
sx q[3];
rz(2.1029162) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6240876) q[0];
sx q[0];
rz(-0.085260304) q[0];
sx q[0];
rz(0.10661539) q[0];
rz(0.34635776) q[1];
sx q[1];
rz(-2.296591) q[1];
sx q[1];
rz(1.9812298) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24460565) q[0];
sx q[0];
rz(-1.7443313) q[0];
sx q[0];
rz(-0.24994295) q[0];
x q[1];
rz(2.3290655) q[2];
sx q[2];
rz(-0.90355325) q[2];
sx q[2];
rz(0.086712547) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.47487709) q[1];
sx q[1];
rz(-1.2704388) q[1];
sx q[1];
rz(-1.9431861) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87507803) q[3];
sx q[3];
rz(-2.4949673) q[3];
sx q[3];
rz(-2.2602606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0059263) q[2];
sx q[2];
rz(-0.29834193) q[2];
sx q[2];
rz(-3.0653817) q[2];
rz(2.5557319) q[3];
sx q[3];
rz(-1.9867089) q[3];
sx q[3];
rz(-1.7990254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0214486) q[0];
sx q[0];
rz(-1.3055389) q[0];
sx q[0];
rz(0.35010499) q[0];
rz(0.95589751) q[1];
sx q[1];
rz(-1.2799542) q[1];
sx q[1];
rz(1.9116481) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0098388) q[0];
sx q[0];
rz(-1.810605) q[0];
sx q[0];
rz(1.1935704) q[0];
rz(-pi) q[1];
rz(-1.3558055) q[2];
sx q[2];
rz(-2.9724398) q[2];
sx q[2];
rz(-0.40119888) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1175244) q[1];
sx q[1];
rz(-1.9000016) q[1];
sx q[1];
rz(0.19964053) q[1];
rz(-2.941972) q[3];
sx q[3];
rz(-2.1484882) q[3];
sx q[3];
rz(0.22838372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.89796394) q[2];
sx q[2];
rz(-2.882759) q[2];
sx q[2];
rz(-1.6492856) q[2];
rz(-2.7367075) q[3];
sx q[3];
rz(-0.76571524) q[3];
sx q[3];
rz(-2.305472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1355857) q[0];
sx q[0];
rz(-1.0772935) q[0];
sx q[0];
rz(-2.5591922) q[0];
rz(0.68663418) q[1];
sx q[1];
rz(-1.4056987) q[1];
sx q[1];
rz(1.883421) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1863886) q[0];
sx q[0];
rz(-1.135728) q[0];
sx q[0];
rz(3.0032936) q[0];
rz(-1.0354543) q[2];
sx q[2];
rz(-1.8796433) q[2];
sx q[2];
rz(0.18986407) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68738378) q[1];
sx q[1];
rz(-1.6404248) q[1];
sx q[1];
rz(1.6268262) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2518796) q[3];
sx q[3];
rz(-2.5237759) q[3];
sx q[3];
rz(-1.8998673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76576343) q[2];
sx q[2];
rz(-1.148369) q[2];
sx q[2];
rz(1.0180391) q[2];
rz(2.8954519) q[3];
sx q[3];
rz(-1.7459511) q[3];
sx q[3];
rz(1.5181946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4356284) q[0];
sx q[0];
rz(-2.3376597) q[0];
sx q[0];
rz(2.5614118) q[0];
rz(0.14353453) q[1];
sx q[1];
rz(-2.6519471) q[1];
sx q[1];
rz(-0.20763436) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.872179) q[0];
sx q[0];
rz(-2.2996443) q[0];
sx q[0];
rz(-2.6711885) q[0];
x q[1];
rz(2.1499499) q[2];
sx q[2];
rz(-1.8657547) q[2];
sx q[2];
rz(2.7911012) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5112111) q[1];
sx q[1];
rz(-0.89511725) q[1];
sx q[1];
rz(0.33501321) q[1];
x q[2];
rz(-2.8945699) q[3];
sx q[3];
rz(-2.7678856) q[3];
sx q[3];
rz(-1.6173687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2988854) q[2];
sx q[2];
rz(-1.8536114) q[2];
sx q[2];
rz(2.5353954) q[2];
rz(-2.4541564) q[3];
sx q[3];
rz(-1.1339374) q[3];
sx q[3];
rz(-2.4351951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6523478) q[0];
sx q[0];
rz(-0.027996538) q[0];
sx q[0];
rz(-2.0835173) q[0];
rz(0.10969133) q[1];
sx q[1];
rz(-1.1166162) q[1];
sx q[1];
rz(-1.441997) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2588651) q[0];
sx q[0];
rz(-1.8804714) q[0];
sx q[0];
rz(0.41559269) q[0];
rz(-1.2953561) q[2];
sx q[2];
rz(-2.0552962) q[2];
sx q[2];
rz(0.4415919) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7322526) q[1];
sx q[1];
rz(-2.2870018) q[1];
sx q[1];
rz(2.2316037) q[1];
rz(-2.9061716) q[3];
sx q[3];
rz(-1.306433) q[3];
sx q[3];
rz(-0.89371577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.66120061) q[2];
sx q[2];
rz(-2.9471687) q[2];
sx q[2];
rz(1.9332168) q[2];
rz(-0.66323534) q[3];
sx q[3];
rz(-1.4570718) q[3];
sx q[3];
rz(-1.9251582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97063589) q[0];
sx q[0];
rz(-2.1747776) q[0];
sx q[0];
rz(3.048625) q[0];
rz(-1.8611106) q[1];
sx q[1];
rz(-2.4130776) q[1];
sx q[1];
rz(3.0063937) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7057892) q[0];
sx q[0];
rz(-3.0717141) q[0];
sx q[0];
rz(2.4326434) q[0];
rz(-pi) q[1];
rz(0.2580169) q[2];
sx q[2];
rz(-1.4224367) q[2];
sx q[2];
rz(-2.6515863) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7635423) q[1];
sx q[1];
rz(-1.5494746) q[1];
sx q[1];
rz(0.2539907) q[1];
x q[2];
rz(-1.1208833) q[3];
sx q[3];
rz(-0.16928798) q[3];
sx q[3];
rz(-1.5811046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4325503) q[2];
sx q[2];
rz(-1.21864) q[2];
sx q[2];
rz(2.3700628) q[2];
rz(-2.6500474) q[3];
sx q[3];
rz(-1.2794269) q[3];
sx q[3];
rz(2.5468723) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5563357) q[0];
sx q[0];
rz(-2.519643) q[0];
sx q[0];
rz(-1.1248032) q[0];
rz(1.8428165) q[1];
sx q[1];
rz(-0.61538428) q[1];
sx q[1];
rz(-2.3769456) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.619304) q[0];
sx q[0];
rz(-1.4656656) q[0];
sx q[0];
rz(-1.6572857) q[0];
rz(-pi) q[1];
rz(0.89057335) q[2];
sx q[2];
rz(-1.9720417) q[2];
sx q[2];
rz(0.97637343) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.54658871) q[1];
sx q[1];
rz(-1.4440795) q[1];
sx q[1];
rz(-1.7751883) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9493136) q[3];
sx q[3];
rz(-2.0928185) q[3];
sx q[3];
rz(-1.6917563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3760959) q[2];
sx q[2];
rz(-1.8509879) q[2];
sx q[2];
rz(0.20720227) q[2];
rz(-2.1616705) q[3];
sx q[3];
rz(-1.1332952) q[3];
sx q[3];
rz(0.19237147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1205263) q[0];
sx q[0];
rz(-1.4561894) q[0];
sx q[0];
rz(-0.86984632) q[0];
rz(0.5008685) q[1];
sx q[1];
rz(-2.905838) q[1];
sx q[1];
rz(0.9314608) q[1];
rz(1.6899213) q[2];
sx q[2];
rz(-1.7075734) q[2];
sx q[2];
rz(1.3592958) q[2];
rz(0.23811447) q[3];
sx q[3];
rz(-1.4128216) q[3];
sx q[3];
rz(-2.9924389) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
