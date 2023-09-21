OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72397435) q[0];
sx q[0];
rz(-1.6516049) q[0];
sx q[0];
rz(0.93044257) q[0];
rz(0.62970495) q[1];
sx q[1];
rz(4.2760744) q[1];
sx q[1];
rz(8.3174336) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9288899) q[0];
sx q[0];
rz(-1.2455997) q[0];
sx q[0];
rz(-2.4341499) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0385752) q[2];
sx q[2];
rz(-0.79471171) q[2];
sx q[2];
rz(2.6796535) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2335637) q[1];
sx q[1];
rz(-1.6828013) q[1];
sx q[1];
rz(1.8506552) q[1];
rz(-pi) q[2];
rz(2.0290306) q[3];
sx q[3];
rz(-2.4535865) q[3];
sx q[3];
rz(-1.0109166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.063623108) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(-1.3280274) q[2];
rz(-0.32087457) q[3];
sx q[3];
rz(-2.1556373) q[3];
sx q[3];
rz(3.0096171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6593453) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(-2.2609718) q[0];
rz(-1.2940787) q[1];
sx q[1];
rz(-0.41795119) q[1];
sx q[1];
rz(-0.81726384) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0016804455) q[0];
sx q[0];
rz(-1.7117371) q[0];
sx q[0];
rz(-1.9190448) q[0];
rz(-pi) q[1];
rz(-1.0802644) q[2];
sx q[2];
rz(-1.763952) q[2];
sx q[2];
rz(-2.1527388) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4546928) q[1];
sx q[1];
rz(-0.94140879) q[1];
sx q[1];
rz(-1.8886186) q[1];
rz(-pi) q[2];
rz(-0.18873429) q[3];
sx q[3];
rz(-1.1047603) q[3];
sx q[3];
rz(-0.97944328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2237504) q[2];
sx q[2];
rz(-0.68085256) q[2];
sx q[2];
rz(-2.7775653) q[2];
rz(-2.1552127) q[3];
sx q[3];
rz(-1.4168408) q[3];
sx q[3];
rz(1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7746975) q[0];
sx q[0];
rz(-0.8323454) q[0];
sx q[0];
rz(-0.96631518) q[0];
rz(2.9486588) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(-1.4470709) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4149949) q[0];
sx q[0];
rz(-1.5532171) q[0];
sx q[0];
rz(0.21501644) q[0];
rz(1.8615396) q[2];
sx q[2];
rz(-1.4915407) q[2];
sx q[2];
rz(2.5153164) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6009532) q[1];
sx q[1];
rz(-0.97266957) q[1];
sx q[1];
rz(-1.1126493) q[1];
x q[2];
rz(-2.0855911) q[3];
sx q[3];
rz(-2.6714973) q[3];
sx q[3];
rz(1.521829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43859279) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(1.1052216) q[2];
rz(-2.3953719) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(2.1658649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0063909) q[0];
sx q[0];
rz(-0.52440301) q[0];
sx q[0];
rz(1.4659457) q[0];
rz(-2.8566467) q[1];
sx q[1];
rz(-1.0703215) q[1];
sx q[1];
rz(-2.7526061) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8143512) q[0];
sx q[0];
rz(-1.8530493) q[0];
sx q[0];
rz(1.809657) q[0];
rz(-pi) q[1];
rz(1.9289891) q[2];
sx q[2];
rz(-2.6274523) q[2];
sx q[2];
rz(-0.51353031) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0491838) q[1];
sx q[1];
rz(-2.1158943) q[1];
sx q[1];
rz(1.070302) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3652521) q[3];
sx q[3];
rz(-0.83344978) q[3];
sx q[3];
rz(1.8635686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(-2.6848865) q[2];
rz(1.5151954) q[3];
sx q[3];
rz(-0.94907343) q[3];
sx q[3];
rz(2.8592498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91963768) q[0];
sx q[0];
rz(-2.0018405) q[0];
sx q[0];
rz(2.7815681) q[0];
rz(0.64741627) q[1];
sx q[1];
rz(-1.5292239) q[1];
sx q[1];
rz(2.6470851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7617103) q[0];
sx q[0];
rz(-2.8344791) q[0];
sx q[0];
rz(-2.7193927) q[0];
x q[1];
rz(3.102166) q[2];
sx q[2];
rz(-1.7061503) q[2];
sx q[2];
rz(-0.53265041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.20930418) q[1];
sx q[1];
rz(-0.63289019) q[1];
sx q[1];
rz(-3.0401405) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5962283) q[3];
sx q[3];
rz(-2.0147707) q[3];
sx q[3];
rz(0.80029315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4313724) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(-2.8430856) q[2];
rz(-2.8295637) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(0.63849866) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1453778) q[0];
sx q[0];
rz(-2.4185116) q[0];
sx q[0];
rz(0.19590713) q[0];
rz(0.021082489) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(-1.235199) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7346749) q[0];
sx q[0];
rz(-1.9039246) q[0];
sx q[0];
rz(-1.8963277) q[0];
x q[1];
rz(-2.7487095) q[2];
sx q[2];
rz(-2.3977931) q[2];
sx q[2];
rz(1.5768029) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0022618) q[1];
sx q[1];
rz(-1.1552703) q[1];
sx q[1];
rz(2.5050487) q[1];
rz(2.9052832) q[3];
sx q[3];
rz(-1.4893388) q[3];
sx q[3];
rz(-2.8636275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49332508) q[2];
sx q[2];
rz(-0.92612925) q[2];
sx q[2];
rz(2.7098999) q[2];
rz(1.7290944) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(0.13599642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39111185) q[0];
sx q[0];
rz(-1.1688122) q[0];
sx q[0];
rz(0.11211638) q[0];
rz(2.926459) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(-1.1134061) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33681413) q[0];
sx q[0];
rz(-1.830173) q[0];
sx q[0];
rz(2.74385) q[0];
rz(-pi) q[1];
rz(0.68219296) q[2];
sx q[2];
rz(-0.39794121) q[2];
sx q[2];
rz(-1.8856018) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.73247319) q[1];
sx q[1];
rz(-1.2670244) q[1];
sx q[1];
rz(-2.4227546) q[1];
x q[2];
rz(-0.21344276) q[3];
sx q[3];
rz(-2.0091972) q[3];
sx q[3];
rz(-1.9813117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1404184) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(0.12602885) q[2];
rz(-1.0472939) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(0.70820156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4814608) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(1.6814167) q[0];
rz(1.2449645) q[1];
sx q[1];
rz(-1.0943202) q[1];
sx q[1];
rz(1.2089027) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7891114) q[0];
sx q[0];
rz(-1.4915823) q[0];
sx q[0];
rz(-1.3080025) q[0];
x q[1];
rz(0.86645855) q[2];
sx q[2];
rz(-0.87477113) q[2];
sx q[2];
rz(2.8528086) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.066684494) q[1];
sx q[1];
rz(-1.2411331) q[1];
sx q[1];
rz(-2.3782905) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2743837) q[3];
sx q[3];
rz(-0.65215462) q[3];
sx q[3];
rz(-1.46331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.37844354) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(-1.8219927) q[2];
rz(-0.59213263) q[3];
sx q[3];
rz(-1.7254646) q[3];
sx q[3];
rz(0.035141703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.9086583) q[0];
sx q[0];
rz(-0.52353752) q[0];
sx q[0];
rz(-1.7804902) q[0];
rz(1.2110442) q[1];
sx q[1];
rz(-2.2369604) q[1];
sx q[1];
rz(-2.7499054) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64205326) q[0];
sx q[0];
rz(-1.5649438) q[0];
sx q[0];
rz(0.43453479) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6892745) q[2];
sx q[2];
rz(-1.3903793) q[2];
sx q[2];
rz(0.42186055) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3561331) q[1];
sx q[1];
rz(-1.2411989) q[1];
sx q[1];
rz(2.3258924) q[1];
rz(-2.6158995) q[3];
sx q[3];
rz(-1.5075397) q[3];
sx q[3];
rz(-2.4478108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6212375) q[2];
sx q[2];
rz(-1.3427799) q[2];
sx q[2];
rz(1.8048145) q[2];
rz(-2.3796066) q[3];
sx q[3];
rz(-2.8218994) q[3];
sx q[3];
rz(2.3412162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8005463) q[0];
sx q[0];
rz(-0.30277345) q[0];
sx q[0];
rz(2.5706932) q[0];
rz(-1.7123429) q[1];
sx q[1];
rz(-2.0776904) q[1];
sx q[1];
rz(0.16194078) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25625944) q[0];
sx q[0];
rz(-0.67978871) q[0];
sx q[0];
rz(1.1023561) q[0];
rz(1.759388) q[2];
sx q[2];
rz(-1.57975) q[2];
sx q[2];
rz(0.21536516) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.098617741) q[1];
sx q[1];
rz(-2.9264755) q[1];
sx q[1];
rz(-2.2153562) q[1];
rz(0.057823618) q[3];
sx q[3];
rz(-2.6177546) q[3];
sx q[3];
rz(-0.1334838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.95742115) q[2];
sx q[2];
rz(-1.1051757) q[2];
sx q[2];
rz(-1.4282248) q[2];
rz(1.1994294) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(1.7709581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7006871) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(1.3700925) q[1];
sx q[1];
rz(-2.1961828) q[1];
sx q[1];
rz(-0.97074769) q[1];
rz(0.8926819) q[2];
sx q[2];
rz(-2.953139) q[2];
sx q[2];
rz(2.1072731) q[2];
rz(0.42967038) q[3];
sx q[3];
rz(-1.3702787) q[3];
sx q[3];
rz(-2.7313781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
