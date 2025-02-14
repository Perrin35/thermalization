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
rz(2.4801369) q[0];
sx q[0];
rz(-2.8602726) q[0];
sx q[0];
rz(1.9982279) q[0];
rz(0.65302628) q[1];
sx q[1];
rz(4.9622494) q[1];
sx q[1];
rz(9.4402037) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2020039) q[0];
sx q[0];
rz(-2.1103854) q[0];
sx q[0];
rz(2.5651776) q[0];
rz(-2.562343) q[2];
sx q[2];
rz(-1.0760032) q[2];
sx q[2];
rz(-0.65828568) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.676486) q[1];
sx q[1];
rz(-2.3576405) q[1];
sx q[1];
rz(-2.1539861) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7687092) q[3];
sx q[3];
rz(-1.0413525) q[3];
sx q[3];
rz(0.44482619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36110863) q[2];
sx q[2];
rz(-3.1321654) q[2];
sx q[2];
rz(-1.1040322) q[2];
rz(0.55936724) q[3];
sx q[3];
rz(-2.1909824) q[3];
sx q[3];
rz(-2.2573788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1169432) q[0];
sx q[0];
rz(-0.58498061) q[0];
sx q[0];
rz(2.2351433) q[0];
rz(0.34140423) q[1];
sx q[1];
rz(-2.2661792) q[1];
sx q[1];
rz(2.7934449) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9885555) q[0];
sx q[0];
rz(-1.2619902) q[0];
sx q[0];
rz(1.8304118) q[0];
rz(-2.1402713) q[2];
sx q[2];
rz(-3.0932326) q[2];
sx q[2];
rz(2.0404301) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7469937) q[1];
sx q[1];
rz(-2.1809757) q[1];
sx q[1];
rz(-3.0501409) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8208037) q[3];
sx q[3];
rz(-0.75859447) q[3];
sx q[3];
rz(1.7649253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1657095) q[2];
sx q[2];
rz(-2.7996863) q[2];
sx q[2];
rz(-0.085414097) q[2];
rz(-0.092770569) q[3];
sx q[3];
rz(-0.86920357) q[3];
sx q[3];
rz(-0.87576491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3669325) q[0];
sx q[0];
rz(-2.7800738) q[0];
sx q[0];
rz(2.7969978) q[0];
rz(-1.5917646) q[1];
sx q[1];
rz(-1.4180309) q[1];
sx q[1];
rz(-1.6835015) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2705459) q[0];
sx q[0];
rz(-0.71840175) q[0];
sx q[0];
rz(0.055553546) q[0];
x q[1];
rz(-1.6750653) q[2];
sx q[2];
rz(-0.51926702) q[2];
sx q[2];
rz(2.7291275) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.96839206) q[1];
sx q[1];
rz(-1.3276982) q[1];
sx q[1];
rz(0.17420423) q[1];
rz(2.167554) q[3];
sx q[3];
rz(-1.1263873) q[3];
sx q[3];
rz(0.35154464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.496326) q[2];
sx q[2];
rz(-0.93924773) q[2];
sx q[2];
rz(1.1305031) q[2];
rz(-0.94275236) q[3];
sx q[3];
rz(-2.0944984) q[3];
sx q[3];
rz(1.9345136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4261674) q[0];
sx q[0];
rz(-1.3454477) q[0];
sx q[0];
rz(-1.6834393) q[0];
rz(-2.093105) q[1];
sx q[1];
rz(-1.4920934) q[1];
sx q[1];
rz(-0.47007158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9662358) q[0];
sx q[0];
rz(-2.5585173) q[0];
sx q[0];
rz(-0.52578853) q[0];
rz(-2.0087035) q[2];
sx q[2];
rz(-2.5688969) q[2];
sx q[2];
rz(2.0349791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4973891) q[1];
sx q[1];
rz(-0.70442373) q[1];
sx q[1];
rz(-0.299244) q[1];
rz(-pi) q[2];
rz(0.87155452) q[3];
sx q[3];
rz(-1.0207173) q[3];
sx q[3];
rz(-0.14000237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5992392) q[2];
sx q[2];
rz(-2.3848644) q[2];
sx q[2];
rz(2.2139464) q[2];
rz(1.8574235) q[3];
sx q[3];
rz(-1.7312867) q[3];
sx q[3];
rz(-0.94902432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014378431) q[0];
sx q[0];
rz(-0.40507409) q[0];
sx q[0];
rz(0.48144427) q[0];
rz(-2.1638347) q[1];
sx q[1];
rz(-2.4974186) q[1];
sx q[1];
rz(-2.0668623) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18431379) q[0];
sx q[0];
rz(-1.1536351) q[0];
sx q[0];
rz(2.8431312) q[0];
rz(-pi) q[1];
rz(2.3911209) q[2];
sx q[2];
rz(-2.7402862) q[2];
sx q[2];
rz(-1.4913781) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9323651) q[1];
sx q[1];
rz(-1.3316939) q[1];
sx q[1];
rz(1.9537611) q[1];
rz(-pi) q[2];
rz(1.2541708) q[3];
sx q[3];
rz(-1.9434557) q[3];
sx q[3];
rz(2.0697442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0398728) q[2];
sx q[2];
rz(-1.1168672) q[2];
sx q[2];
rz(-2.5679585) q[2];
rz(-1.6576069) q[3];
sx q[3];
rz(-2.0935757) q[3];
sx q[3];
rz(0.60605961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1178591) q[0];
sx q[0];
rz(-0.25256279) q[0];
sx q[0];
rz(-0.31164393) q[0];
rz(-0.32870865) q[1];
sx q[1];
rz(-1.3361822) q[1];
sx q[1];
rz(-2.4005344) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58068507) q[0];
sx q[0];
rz(-2.0253165) q[0];
sx q[0];
rz(0.88991388) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0799505) q[2];
sx q[2];
rz(-1.3423402) q[2];
sx q[2];
rz(2.1232587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0892798) q[1];
sx q[1];
rz(-1.5075753) q[1];
sx q[1];
rz(2.2526645) q[1];
rz(-pi) q[2];
rz(-0.00087329207) q[3];
sx q[3];
rz(-2.2016618) q[3];
sx q[3];
rz(-1.445595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3662423) q[2];
sx q[2];
rz(-2.5516208) q[2];
sx q[2];
rz(-0.74014202) q[2];
rz(2.7548693) q[3];
sx q[3];
rz(-2.4155278) q[3];
sx q[3];
rz(-0.47851419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9750403) q[0];
sx q[0];
rz(-2.1599025) q[0];
sx q[0];
rz(0.24328406) q[0];
rz(-0.89726204) q[1];
sx q[1];
rz(-1.5829395) q[1];
sx q[1];
rz(0.74403393) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99951852) q[0];
sx q[0];
rz(-1.4878977) q[0];
sx q[0];
rz(-2.7510452) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7626708) q[2];
sx q[2];
rz(-2.1728467) q[2];
sx q[2];
rz(-3.0343891) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8764827) q[1];
sx q[1];
rz(-1.0534673) q[1];
sx q[1];
rz(-0.4267652) q[1];
rz(2.9570368) q[3];
sx q[3];
rz(-2.089644) q[3];
sx q[3];
rz(-0.59901064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1870785) q[2];
sx q[2];
rz(-2.3176471) q[2];
sx q[2];
rz(0.0079060923) q[2];
rz(1.6451969) q[3];
sx q[3];
rz(-3.0683066) q[3];
sx q[3];
rz(0.50905281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024648333) q[0];
sx q[0];
rz(-0.032289676) q[0];
sx q[0];
rz(-3.0049128) q[0];
rz(-1.7928803) q[1];
sx q[1];
rz(-1.2929448) q[1];
sx q[1];
rz(0.50833702) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6240468) q[0];
sx q[0];
rz(-1.4410748) q[0];
sx q[0];
rz(1.172439) q[0];
x q[1];
rz(2.2466564) q[2];
sx q[2];
rz(-1.3351262) q[2];
sx q[2];
rz(-1.6431944) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.52259787) q[1];
sx q[1];
rz(-2.3648089) q[1];
sx q[1];
rz(2.0911699) q[1];
rz(-pi) q[2];
x q[2];
rz(3.062183) q[3];
sx q[3];
rz(-1.4520922) q[3];
sx q[3];
rz(2.1044923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4623744) q[2];
sx q[2];
rz(-2.4269673) q[2];
sx q[2];
rz(0.057961658) q[2];
rz(-0.95311779) q[3];
sx q[3];
rz(-2.0942196) q[3];
sx q[3];
rz(0.63547772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5304831) q[0];
sx q[0];
rz(-1.4169175) q[0];
sx q[0];
rz(-0.86607754) q[0];
rz(-1.3653612) q[1];
sx q[1];
rz(-1.583464) q[1];
sx q[1];
rz(0.80397111) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3063076) q[0];
sx q[0];
rz(-0.87125766) q[0];
sx q[0];
rz(0.37658306) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3525904) q[2];
sx q[2];
rz(-1.2478634) q[2];
sx q[2];
rz(-2.5247521) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6996346) q[1];
sx q[1];
rz(-0.82773877) q[1];
sx q[1];
rz(-0.98299594) q[1];
rz(1.1671806) q[3];
sx q[3];
rz(-1.5909373) q[3];
sx q[3];
rz(2.8668364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.08635252) q[2];
sx q[2];
rz(-2.99282) q[2];
sx q[2];
rz(-1.0521592) q[2];
rz(-0.85938984) q[3];
sx q[3];
rz(-2.342577) q[3];
sx q[3];
rz(-0.94256443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6732442) q[0];
sx q[0];
rz(-2.4788661) q[0];
sx q[0];
rz(-0.41917875) q[0];
rz(-3.1255417) q[1];
sx q[1];
rz(-1.5546067) q[1];
sx q[1];
rz(0.12414653) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9036839) q[0];
sx q[0];
rz(-0.57460898) q[0];
sx q[0];
rz(1.2875071) q[0];
x q[1];
rz(1.8252462) q[2];
sx q[2];
rz(-1.7702603) q[2];
sx q[2];
rz(0.9141038) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.032822306) q[1];
sx q[1];
rz(-1.747393) q[1];
sx q[1];
rz(-2.887421) q[1];
rz(0.92049034) q[3];
sx q[3];
rz(-1.4249885) q[3];
sx q[3];
rz(2.15964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9496256) q[2];
sx q[2];
rz(-2.4032335) q[2];
sx q[2];
rz(2.9730566) q[2];
rz(-1.4847697) q[3];
sx q[3];
rz(-0.46983457) q[3];
sx q[3];
rz(2.7114939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1930502) q[0];
sx q[0];
rz(-0.47946231) q[0];
sx q[0];
rz(-0.23647501) q[0];
rz(-0.82462689) q[1];
sx q[1];
rz(-1.5390479) q[1];
sx q[1];
rz(1.8631757) q[1];
rz(-2.5599418) q[2];
sx q[2];
rz(-0.86514513) q[2];
sx q[2];
rz(-1.9164597) q[2];
rz(3.0807224) q[3];
sx q[3];
rz(-1.0247598) q[3];
sx q[3];
rz(1.8714874) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
