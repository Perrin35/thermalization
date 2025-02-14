OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.4247704) q[0];
sx q[0];
rz(4.7228887) q[0];
sx q[0];
rz(6.9520998) q[0];
rz(-2.7954697) q[1];
sx q[1];
rz(-0.82155138) q[1];
sx q[1];
rz(0.93924826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6990841) q[0];
sx q[0];
rz(-2.6538349) q[0];
sx q[0];
rz(0.91983025) q[0];
x q[1];
rz(3.0947565) q[2];
sx q[2];
rz(-1.9745262) q[2];
sx q[2];
rz(0.22778928) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2722764) q[1];
sx q[1];
rz(-2.6679435) q[1];
sx q[1];
rz(0.95817153) q[1];
rz(-pi) q[2];
rz(-1.2143308) q[3];
sx q[3];
rz(-2.2495396) q[3];
sx q[3];
rz(-1.9744622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5758489) q[2];
sx q[2];
rz(-1.0172903) q[2];
sx q[2];
rz(3.1080833) q[2];
rz(-1.2014028) q[3];
sx q[3];
rz(-1.3591432) q[3];
sx q[3];
rz(1.0777773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1101542) q[0];
sx q[0];
rz(-0.25122508) q[0];
sx q[0];
rz(0.95397368) q[0];
rz(-1.4454449) q[1];
sx q[1];
rz(-2.0868389) q[1];
sx q[1];
rz(1.1711858) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508043) q[0];
sx q[0];
rz(-0.91676869) q[0];
sx q[0];
rz(-0.79933856) q[0];
rz(-pi) q[1];
rz(1.2759802) q[2];
sx q[2];
rz(-0.95446842) q[2];
sx q[2];
rz(2.4957239) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7668768) q[1];
sx q[1];
rz(-1.7991814) q[1];
sx q[1];
rz(1.5414052) q[1];
rz(-pi) q[2];
rz(2.0061139) q[3];
sx q[3];
rz(-1.2655228) q[3];
sx q[3];
rz(-1.0686309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.235432) q[2];
sx q[2];
rz(-1.7110598) q[2];
sx q[2];
rz(1.1492427) q[2];
rz(-3.1084642) q[3];
sx q[3];
rz(-1.5002316) q[3];
sx q[3];
rz(-0.59605014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(1.0572492) q[0];
sx q[0];
rz(-0.79279041) q[0];
sx q[0];
rz(-1.0302011) q[0];
rz(1.239981) q[1];
sx q[1];
rz(-1.3571309) q[1];
sx q[1];
rz(0.99303594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030331197) q[0];
sx q[0];
rz(-1.5227888) q[0];
sx q[0];
rz(2.6859849) q[0];
rz(-0.15056653) q[2];
sx q[2];
rz(-0.67306256) q[2];
sx q[2];
rz(1.7889495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5883623) q[1];
sx q[1];
rz(-2.2892286) q[1];
sx q[1];
rz(-2.3521982) q[1];
rz(-pi) q[2];
rz(-2.7168324) q[3];
sx q[3];
rz(-2.3094607) q[3];
sx q[3];
rz(0.55603851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21732907) q[2];
sx q[2];
rz(-2.0254878) q[2];
sx q[2];
rz(0.35476157) q[2];
rz(-0.99700704) q[3];
sx q[3];
rz(-1.487178) q[3];
sx q[3];
rz(2.2009489) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16911258) q[0];
sx q[0];
rz(-1.4262154) q[0];
sx q[0];
rz(-2.0347563) q[0];
rz(2.2726982) q[1];
sx q[1];
rz(-2.6241701) q[1];
sx q[1];
rz(-0.69721627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1887564) q[0];
sx q[0];
rz(-1.3719014) q[0];
sx q[0];
rz(-2.5752978) q[0];
rz(-0.12733404) q[2];
sx q[2];
rz(-2.4155354) q[2];
sx q[2];
rz(0.51148326) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6706411) q[1];
sx q[1];
rz(-1.9578251) q[1];
sx q[1];
rz(0.3989889) q[1];
x q[2];
rz(-0.02454464) q[3];
sx q[3];
rz(-1.595579) q[3];
sx q[3];
rz(-1.9056232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8640459) q[2];
sx q[2];
rz(-2.9329381) q[2];
sx q[2];
rz(0.29402688) q[2];
rz(-1.0634408) q[3];
sx q[3];
rz(-1.4592417) q[3];
sx q[3];
rz(2.3991876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4251637) q[0];
sx q[0];
rz(-0.18121457) q[0];
sx q[0];
rz(2.678405) q[0];
rz(-0.40395346) q[1];
sx q[1];
rz(-1.5084167) q[1];
sx q[1];
rz(0.24872669) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8578019) q[0];
sx q[0];
rz(-1.3732855) q[0];
sx q[0];
rz(3.1398612) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.727365) q[2];
sx q[2];
rz(-0.58097208) q[2];
sx q[2];
rz(1.5805949) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.675093) q[1];
sx q[1];
rz(-0.50216253) q[1];
sx q[1];
rz(-2.2999022) q[1];
x q[2];
rz(-2.4047923) q[3];
sx q[3];
rz(-2.1016788) q[3];
sx q[3];
rz(-0.30217043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2948461) q[2];
sx q[2];
rz(-1.8241901) q[2];
sx q[2];
rz(1.741629) q[2];
rz(1.8585662) q[3];
sx q[3];
rz(-1.3834407) q[3];
sx q[3];
rz(-2.9156901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.95752174) q[0];
sx q[0];
rz(-2.3029843) q[0];
sx q[0];
rz(2.7958909) q[0];
rz(0.050994571) q[1];
sx q[1];
rz(-1.9146999) q[1];
sx q[1];
rz(-2.3695703) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2564023) q[0];
sx q[0];
rz(-1.2202383) q[0];
sx q[0];
rz(-1.3023443) q[0];
x q[1];
rz(0.096944158) q[2];
sx q[2];
rz(-2.2205995) q[2];
sx q[2];
rz(0.61496269) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.0061969697) q[1];
sx q[1];
rz(-0.67402285) q[1];
sx q[1];
rz(-2.1042906) q[1];
rz(-pi) q[2];
rz(2.0421626) q[3];
sx q[3];
rz(-2.5623218) q[3];
sx q[3];
rz(1.8806277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2766075) q[2];
sx q[2];
rz(-1.425068) q[2];
sx q[2];
rz(2.9841606) q[2];
rz(-0.43846798) q[3];
sx q[3];
rz(-2.4003568) q[3];
sx q[3];
rz(0.95611519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0601783) q[0];
sx q[0];
rz(-2.0745451) q[0];
sx q[0];
rz(-9/(11*pi)) q[0];
rz(-2.5632437) q[1];
sx q[1];
rz(-1.3713505) q[1];
sx q[1];
rz(-0.15301212) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3653202) q[0];
sx q[0];
rz(-2.0663102) q[0];
sx q[0];
rz(2.1644664) q[0];
rz(0.084916755) q[2];
sx q[2];
rz(-1.6667637) q[2];
sx q[2];
rz(-1.7964448) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.12100231) q[1];
sx q[1];
rz(-0.7210702) q[1];
sx q[1];
rz(-2.4929558) q[1];
rz(-0.97665321) q[3];
sx q[3];
rz(-0.53721957) q[3];
sx q[3];
rz(1.9013167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.962062) q[2];
sx q[2];
rz(-0.099055812) q[2];
sx q[2];
rz(2.8774101) q[2];
rz(-1.3029441) q[3];
sx q[3];
rz(-1.3774201) q[3];
sx q[3];
rz(2.0343659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1846979) q[0];
sx q[0];
rz(-0.95320025) q[0];
sx q[0];
rz(1.6860929) q[0];
rz(0.9616583) q[1];
sx q[1];
rz(-1.3788297) q[1];
sx q[1];
rz(-2.2419194) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6754705) q[0];
sx q[0];
rz(-0.5044901) q[0];
sx q[0];
rz(-2.6928508) q[0];
rz(-pi) q[1];
rz(-1.8223962) q[2];
sx q[2];
rz(-2.4008022) q[2];
sx q[2];
rz(-1.0902001) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4845061) q[1];
sx q[1];
rz(-2.5528347) q[1];
sx q[1];
rz(1.4935054) q[1];
rz(-1.353823) q[3];
sx q[3];
rz(-1.5343175) q[3];
sx q[3];
rz(-2.77751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6962894) q[2];
sx q[2];
rz(-0.60089198) q[2];
sx q[2];
rz(0.75330934) q[2];
rz(0.2937915) q[3];
sx q[3];
rz(-2.0132422) q[3];
sx q[3];
rz(1.1118579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488572) q[0];
sx q[0];
rz(-0.98965544) q[0];
sx q[0];
rz(-2.2897172) q[0];
rz(1.4082255) q[1];
sx q[1];
rz(-1.4829166) q[1];
sx q[1];
rz(2.7588989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4538759) q[0];
sx q[0];
rz(-1.9864489) q[0];
sx q[0];
rz(1.835653) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23613249) q[2];
sx q[2];
rz(-0.78401596) q[2];
sx q[2];
rz(1.2003984) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.574461) q[1];
sx q[1];
rz(-1.2700915) q[1];
sx q[1];
rz(1.6733132) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15573536) q[3];
sx q[3];
rz(-2.6504575) q[3];
sx q[3];
rz(0.30065003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7644299) q[2];
sx q[2];
rz(-2.1775776) q[2];
sx q[2];
rz(2.5386179) q[2];
rz(2.073334) q[3];
sx q[3];
rz(-1.7030741) q[3];
sx q[3];
rz(-2.300613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0076440796) q[0];
sx q[0];
rz(-1.4643865) q[0];
sx q[0];
rz(0.35368791) q[0];
rz(2.0091281) q[1];
sx q[1];
rz(-0.9681038) q[1];
sx q[1];
rz(-0.38945928) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4492466) q[0];
sx q[0];
rz(-2.5857537) q[0];
sx q[0];
rz(0.44958322) q[0];
rz(-pi) q[1];
rz(1.816941) q[2];
sx q[2];
rz(-1.4708286) q[2];
sx q[2];
rz(2.113518) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88122565) q[1];
sx q[1];
rz(-1.1133476) q[1];
sx q[1];
rz(1.264099) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0515303) q[3];
sx q[3];
rz(-0.52999338) q[3];
sx q[3];
rz(1.8217721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0103717) q[2];
sx q[2];
rz(-1.5067357) q[2];
sx q[2];
rz(-1.1466675) q[2];
rz(1.6049339) q[3];
sx q[3];
rz(-2.6585572) q[3];
sx q[3];
rz(0.35227942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1491886) q[0];
sx q[0];
rz(-2.438899) q[0];
sx q[0];
rz(-2.6371523) q[0];
rz(-0.06012499) q[1];
sx q[1];
rz(-0.45777121) q[1];
sx q[1];
rz(-0.4578185) q[1];
rz(-1.6788775) q[2];
sx q[2];
rz(-1.6440132) q[2];
sx q[2];
rz(-3.0655412) q[2];
rz(-3.1067228) q[3];
sx q[3];
rz(-2.2449319) q[3];
sx q[3];
rz(1.9064376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
