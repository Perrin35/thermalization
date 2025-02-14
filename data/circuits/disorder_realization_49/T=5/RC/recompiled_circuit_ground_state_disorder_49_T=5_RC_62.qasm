OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9743118) q[0];
sx q[0];
rz(-0.33742961) q[0];
sx q[0];
rz(-2.9396102) q[0];
rz(2.1972411) q[1];
sx q[1];
rz(-0.75466067) q[1];
sx q[1];
rz(-1.6973629) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026469507) q[0];
sx q[0];
rz(-1.7805003) q[0];
sx q[0];
rz(-0.50499495) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0336122) q[2];
sx q[2];
rz(-1.3657346) q[2];
sx q[2];
rz(-0.16879665) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7601499) q[1];
sx q[1];
rz(-2.0408148) q[1];
sx q[1];
rz(2.2862741) q[1];
rz(-pi) q[2];
rz(-1.6891278) q[3];
sx q[3];
rz(-0.36709309) q[3];
sx q[3];
rz(-1.5637507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.28295383) q[2];
sx q[2];
rz(-3.0573461) q[2];
sx q[2];
rz(-1.6072744) q[2];
rz(0.18680799) q[3];
sx q[3];
rz(-0.45972937) q[3];
sx q[3];
rz(1.648858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9635791) q[0];
sx q[0];
rz(-2.3602965) q[0];
sx q[0];
rz(0.70469967) q[0];
rz(2.1251382) q[1];
sx q[1];
rz(-1.1926788) q[1];
sx q[1];
rz(1.1526398) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3160313) q[0];
sx q[0];
rz(-1.7411147) q[0];
sx q[0];
rz(-1.0206778) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5637758) q[2];
sx q[2];
rz(-1.0989099) q[2];
sx q[2];
rz(-3.0834215) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1543442) q[1];
sx q[1];
rz(-1.0030522) q[1];
sx q[1];
rz(0.61056925) q[1];
rz(0.14039881) q[3];
sx q[3];
rz(-0.75274668) q[3];
sx q[3];
rz(-1.6209768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6013201) q[2];
sx q[2];
rz(-0.86869621) q[2];
sx q[2];
rz(-1.8216088) q[2];
rz(-3.0705304) q[3];
sx q[3];
rz(-3.0607405) q[3];
sx q[3];
rz(-2.0474056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0367301) q[0];
sx q[0];
rz(-2.7920089) q[0];
sx q[0];
rz(0.87376755) q[0];
rz(1.8050487) q[1];
sx q[1];
rz(-0.68160325) q[1];
sx q[1];
rz(-0.85493404) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.076043) q[0];
sx q[0];
rz(-0.67971715) q[0];
sx q[0];
rz(-1.2290086) q[0];
rz(-pi) q[1];
rz(0.31051113) q[2];
sx q[2];
rz(-1.1689671) q[2];
sx q[2];
rz(2.7257811) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1347772) q[1];
sx q[1];
rz(-1.2689277) q[1];
sx q[1];
rz(1.4567489) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35936648) q[3];
sx q[3];
rz(-2.1443261) q[3];
sx q[3];
rz(1.9952967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9285589) q[2];
sx q[2];
rz(-2.010767) q[2];
sx q[2];
rz(2.0862759) q[2];
rz(2.1999551) q[3];
sx q[3];
rz(-2.0703273) q[3];
sx q[3];
rz(0.94282237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.68614352) q[0];
sx q[0];
rz(-0.84857714) q[0];
sx q[0];
rz(-0.69547478) q[0];
rz(1.4675568) q[1];
sx q[1];
rz(-1.8753884) q[1];
sx q[1];
rz(-0.094559018) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.031700121) q[0];
sx q[0];
rz(-1.6188834) q[0];
sx q[0];
rz(1.2926433) q[0];
x q[1];
rz(-2.2945061) q[2];
sx q[2];
rz(-1.0930335) q[2];
sx q[2];
rz(1.0003045) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70909158) q[1];
sx q[1];
rz(-1.3729665) q[1];
sx q[1];
rz(1.2835931) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2572722) q[3];
sx q[3];
rz(-2.4289968) q[3];
sx q[3];
rz(2.7698776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62119421) q[2];
sx q[2];
rz(-0.79402557) q[2];
sx q[2];
rz(0.76148477) q[2];
rz(2.9705808) q[3];
sx q[3];
rz(-1.8595502) q[3];
sx q[3];
rz(0.096693501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3068202) q[0];
sx q[0];
rz(-0.67245317) q[0];
sx q[0];
rz(-2.6493454) q[0];
rz(-1.3193839) q[1];
sx q[1];
rz(-1.4232114) q[1];
sx q[1];
rz(0.48670235) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20535417) q[0];
sx q[0];
rz(-1.9689318) q[0];
sx q[0];
rz(-1.7004844) q[0];
rz(-pi) q[1];
rz(-2.3721061) q[2];
sx q[2];
rz(-1.2956007) q[2];
sx q[2];
rz(1.9857032) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4449205) q[1];
sx q[1];
rz(-1.2081523) q[1];
sx q[1];
rz(-1.1128694) q[1];
rz(-2.0858889) q[3];
sx q[3];
rz(-1.7755701) q[3];
sx q[3];
rz(1.0398454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95185602) q[2];
sx q[2];
rz(-2.2978013) q[2];
sx q[2];
rz(1.0978511) q[2];
rz(2.3073933) q[3];
sx q[3];
rz(-2.4662377) q[3];
sx q[3];
rz(-1.437291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6247691) q[0];
sx q[0];
rz(-2.5747445) q[0];
sx q[0];
rz(-0.045850642) q[0];
rz(2.397873) q[1];
sx q[1];
rz(-0.36879483) q[1];
sx q[1];
rz(3.0812982) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1865538) q[0];
sx q[0];
rz(-0.56772029) q[0];
sx q[0];
rz(-0.7132775) q[0];
rz(-1.1072537) q[2];
sx q[2];
rz(-2.3899979) q[2];
sx q[2];
rz(-1.3857402) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0963147) q[1];
sx q[1];
rz(-1.244925) q[1];
sx q[1];
rz(0.23479825) q[1];
rz(-pi) q[2];
rz(-1.6894108) q[3];
sx q[3];
rz(-2.4595692) q[3];
sx q[3];
rz(-0.31957808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9149949) q[2];
sx q[2];
rz(-1.3264341) q[2];
sx q[2];
rz(0.9787406) q[2];
rz(-1.2616875) q[3];
sx q[3];
rz(-1.0190957) q[3];
sx q[3];
rz(0.88920465) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9642445) q[0];
sx q[0];
rz(-1.526399) q[0];
sx q[0];
rz(-0.88441315) q[0];
rz(-0.85909596) q[1];
sx q[1];
rz(-2.0180549) q[1];
sx q[1];
rz(-2.9335847) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7743083) q[0];
sx q[0];
rz(-1.9242263) q[0];
sx q[0];
rz(-2.6103781) q[0];
rz(-2.5701373) q[2];
sx q[2];
rz(-0.68585912) q[2];
sx q[2];
rz(1.4446007) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1237202) q[1];
sx q[1];
rz(-1.3062638) q[1];
sx q[1];
rz(0.912028) q[1];
rz(-1.0674123) q[3];
sx q[3];
rz(-1.5725822) q[3];
sx q[3];
rz(-1.3942277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7615243) q[2];
sx q[2];
rz(-2.2218349) q[2];
sx q[2];
rz(-2.5999787) q[2];
rz(1.2540865) q[3];
sx q[3];
rz(-1.5518291) q[3];
sx q[3];
rz(-0.90887535) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8054473) q[0];
sx q[0];
rz(-2.8271524) q[0];
sx q[0];
rz(-1.0910777) q[0];
rz(0.82487851) q[1];
sx q[1];
rz(-0.42366091) q[1];
sx q[1];
rz(-1.0728015) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43050736) q[0];
sx q[0];
rz(-0.46784624) q[0];
sx q[0];
rz(-2.7070295) q[0];
rz(1.593441) q[2];
sx q[2];
rz(-1.5839108) q[2];
sx q[2];
rz(0.85044059) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78847105) q[1];
sx q[1];
rz(-1.7040729) q[1];
sx q[1];
rz(2.8972096) q[1];
x q[2];
rz(1.1680702) q[3];
sx q[3];
rz(-1.793705) q[3];
sx q[3];
rz(-2.9532578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9020646) q[2];
sx q[2];
rz(-2.5967279) q[2];
sx q[2];
rz(5/(6*pi)) q[2];
rz(-0.15051633) q[3];
sx q[3];
rz(-2.033332) q[3];
sx q[3];
rz(0.36293852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4758258) q[0];
sx q[0];
rz(-0.095223991) q[0];
sx q[0];
rz(1.6640523) q[0];
rz(0.80329576) q[1];
sx q[1];
rz(-1.3731615) q[1];
sx q[1];
rz(-2.1243748) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4519891) q[0];
sx q[0];
rz(-1.7740806) q[0];
sx q[0];
rz(-3.0102121) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46131409) q[2];
sx q[2];
rz(-1.3405356) q[2];
sx q[2];
rz(2.3264937) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0086144) q[1];
sx q[1];
rz(-0.45615754) q[1];
sx q[1];
rz(-2.4309979) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8152596) q[3];
sx q[3];
rz(-1.8584972) q[3];
sx q[3];
rz(0.74244754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9251755) q[2];
sx q[2];
rz(-1.5745796) q[2];
sx q[2];
rz(2.4851921) q[2];
rz(2.9260855) q[3];
sx q[3];
rz(-1.7301205) q[3];
sx q[3];
rz(-1.6925252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0553174) q[0];
sx q[0];
rz(-1.3441939) q[0];
sx q[0];
rz(-2.2917746) q[0];
rz(-0.96018106) q[1];
sx q[1];
rz(-0.4147059) q[1];
sx q[1];
rz(0.93920952) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612027) q[0];
sx q[0];
rz(-1.3074669) q[0];
sx q[0];
rz(-2.3376746) q[0];
x q[1];
rz(-1.9853815) q[2];
sx q[2];
rz(-1.5778981) q[2];
sx q[2];
rz(-2.9455001) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.997949) q[1];
sx q[1];
rz(-0.74926361) q[1];
sx q[1];
rz(0.56203385) q[1];
x q[2];
rz(-2.0980832) q[3];
sx q[3];
rz(-1.1498972) q[3];
sx q[3];
rz(-0.43779201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5822997) q[2];
sx q[2];
rz(-1.7377995) q[2];
sx q[2];
rz(2.5913008) q[2];
rz(-0.127921) q[3];
sx q[3];
rz(-3.0030799) q[3];
sx q[3];
rz(-2.1730455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6561103) q[0];
sx q[0];
rz(-0.67018296) q[0];
sx q[0];
rz(-0.67076587) q[0];
rz(0.35519629) q[1];
sx q[1];
rz(-1.6658446) q[1];
sx q[1];
rz(2.7459941) q[1];
rz(1.4305461) q[2];
sx q[2];
rz(-1.7127019) q[2];
sx q[2];
rz(0.4455926) q[2];
rz(1.1200503) q[3];
sx q[3];
rz(-1.7786296) q[3];
sx q[3];
rz(2.8854388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
