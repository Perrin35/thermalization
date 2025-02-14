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
rz(0.20198241) q[0];
rz(2.1972411) q[1];
sx q[1];
rz(-0.75466067) q[1];
sx q[1];
rz(1.4442297) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2371191) q[0];
sx q[0];
rz(-2.5982862) q[0];
sx q[0];
rz(-0.41443698) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0927222) q[2];
sx q[2];
rz(-2.9101924) q[2];
sx q[2];
rz(-0.3203985) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7601499) q[1];
sx q[1];
rz(-1.1007778) q[1];
sx q[1];
rz(-0.85531855) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9355447) q[3];
sx q[3];
rz(-1.6131796) q[3];
sx q[3];
rz(-0.11755951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28295383) q[2];
sx q[2];
rz(-3.0573461) q[2];
sx q[2];
rz(1.6072744) q[2];
rz(-0.18680799) q[3];
sx q[3];
rz(-0.45972937) q[3];
sx q[3];
rz(-1.648858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1780136) q[0];
sx q[0];
rz(-2.3602965) q[0];
sx q[0];
rz(-2.436893) q[0];
rz(-2.1251382) q[1];
sx q[1];
rz(-1.1926788) q[1];
sx q[1];
rz(1.9889529) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3160313) q[0];
sx q[0];
rz(-1.7411147) q[0];
sx q[0];
rz(-2.1209149) q[0];
rz(2.5637758) q[2];
sx q[2];
rz(-2.0426828) q[2];
sx q[2];
rz(3.0834215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1543442) q[1];
sx q[1];
rz(-2.1385405) q[1];
sx q[1];
rz(0.61056925) q[1];
x q[2];
rz(-0.14039881) q[3];
sx q[3];
rz(-2.388846) q[3];
sx q[3];
rz(1.5206159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54027259) q[2];
sx q[2];
rz(-2.2728964) q[2];
sx q[2];
rz(1.8216088) q[2];
rz(0.071062239) q[3];
sx q[3];
rz(-0.080852121) q[3];
sx q[3];
rz(2.0474056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1048626) q[0];
sx q[0];
rz(-0.34958378) q[0];
sx q[0];
rz(0.87376755) q[0];
rz(-1.8050487) q[1];
sx q[1];
rz(-0.68160325) q[1];
sx q[1];
rz(0.85493404) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4945472) q[0];
sx q[0];
rz(-0.93699199) q[0];
sx q[0];
rz(0.26453544) q[0];
rz(-2.8310815) q[2];
sx q[2];
rz(-1.1689671) q[2];
sx q[2];
rz(2.7257811) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3745649) q[1];
sx q[1];
rz(-0.32207707) q[1];
sx q[1];
rz(-0.35038553) q[1];
x q[2];
rz(1.072238) q[3];
sx q[3];
rz(-0.66594687) q[3];
sx q[3];
rz(1.3897105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2130337) q[2];
sx q[2];
rz(-1.1308257) q[2];
sx q[2];
rz(-2.0862759) q[2];
rz(2.1999551) q[3];
sx q[3];
rz(-2.0703273) q[3];
sx q[3];
rz(-2.1987703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4554491) q[0];
sx q[0];
rz(-2.2930155) q[0];
sx q[0];
rz(-0.69547478) q[0];
rz(-1.6740359) q[1];
sx q[1];
rz(-1.2662042) q[1];
sx q[1];
rz(-3.0470336) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1098925) q[0];
sx q[0];
rz(-1.5227093) q[0];
sx q[0];
rz(-1.8489494) q[0];
rz(-pi) q[1];
rz(-0.60463011) q[2];
sx q[2];
rz(-2.1994091) q[2];
sx q[2];
rz(-0.18456799) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8670204) q[1];
sx q[1];
rz(-2.7944075) q[1];
sx q[1];
rz(2.1866261) q[1];
rz(0.26043268) q[3];
sx q[3];
rz(-0.89959252) q[3];
sx q[3];
rz(0.77653054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5203984) q[2];
sx q[2];
rz(-2.3475671) q[2];
sx q[2];
rz(2.3801079) q[2];
rz(-2.9705808) q[3];
sx q[3];
rz(-1.2820425) q[3];
sx q[3];
rz(0.096693501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8347725) q[0];
sx q[0];
rz(-0.67245317) q[0];
sx q[0];
rz(-2.6493454) q[0];
rz(-1.8222088) q[1];
sx q[1];
rz(-1.4232114) q[1];
sx q[1];
rz(-0.48670235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20535417) q[0];
sx q[0];
rz(-1.1726609) q[0];
sx q[0];
rz(-1.7004844) q[0];
x q[1];
rz(-0.76948659) q[2];
sx q[2];
rz(-1.2956007) q[2];
sx q[2];
rz(1.1558895) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0943739) q[1];
sx q[1];
rz(-1.9969121) q[1];
sx q[1];
rz(2.7414118) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23426849) q[3];
sx q[3];
rz(-2.0741012) q[3];
sx q[3];
rz(0.64555321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1897366) q[2];
sx q[2];
rz(-0.8437914) q[2];
sx q[2];
rz(-2.0437415) q[2];
rz(-2.3073933) q[3];
sx q[3];
rz(-0.67535496) q[3];
sx q[3];
rz(-1.437291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6247691) q[0];
sx q[0];
rz(-2.5747445) q[0];
sx q[0];
rz(0.045850642) q[0];
rz(-2.397873) q[1];
sx q[1];
rz(-2.7727978) q[1];
sx q[1];
rz(-0.060294453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1865538) q[0];
sx q[0];
rz(-0.56772029) q[0];
sx q[0];
rz(0.7132775) q[0];
rz(-0.39581516) q[2];
sx q[2];
rz(-2.2278053) q[2];
sx q[2];
rz(-0.78570494) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44908811) q[1];
sx q[1];
rz(-1.793024) q[1];
sx q[1];
rz(1.9051993) q[1];
rz(-2.2493717) q[3];
sx q[3];
rz(-1.4961317) q[3];
sx q[3];
rz(-1.1589662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.22659773) q[2];
sx q[2];
rz(-1.8151585) q[2];
sx q[2];
rz(0.9787406) q[2];
rz(-1.8799051) q[3];
sx q[3];
rz(-1.0190957) q[3];
sx q[3];
rz(2.252388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17734811) q[0];
sx q[0];
rz(-1.6151936) q[0];
sx q[0];
rz(2.2571795) q[0];
rz(2.2824967) q[1];
sx q[1];
rz(-2.0180549) q[1];
sx q[1];
rz(0.20800796) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0028864127) q[0];
sx q[0];
rz(-2.0660668) q[0];
sx q[0];
rz(1.1664835) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9875097) q[2];
sx q[2];
rz(-1.0089968) q[2];
sx q[2];
rz(2.3903008) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.87863786) q[1];
sx q[1];
rz(-2.4391101) q[1];
sx q[1];
rz(1.9874057) q[1];
rz(0.0020387928) q[3];
sx q[3];
rz(-1.0674132) q[3];
sx q[3];
rz(0.17755213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38006833) q[2];
sx q[2];
rz(-0.91975776) q[2];
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
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8054473) q[0];
sx q[0];
rz(-0.31444028) q[0];
sx q[0];
rz(-1.0910777) q[0];
rz(0.82487851) q[1];
sx q[1];
rz(-2.7179317) q[1];
sx q[1];
rz(1.0728015) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7110853) q[0];
sx q[0];
rz(-2.6737464) q[0];
sx q[0];
rz(-2.7070295) q[0];
rz(-0.013117803) q[2];
sx q[2];
rz(-1.5934391) q[2];
sx q[2];
rz(-0.72005872) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.271928) q[1];
sx q[1];
rz(-2.8638693) q[1];
sx q[1];
rz(-2.6355993) q[1];
rz(-pi) q[2];
rz(-0.24157584) q[3];
sx q[3];
rz(-1.178587) q[3];
sx q[3];
rz(1.6652312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.23952809) q[2];
sx q[2];
rz(-0.54486474) q[2];
sx q[2];
rz(-2.876335) q[2];
rz(-0.15051633) q[3];
sx q[3];
rz(-2.033332) q[3];
sx q[3];
rz(-2.7786541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85452138) q[0];
sx q[0];
rz(-1.6994564) q[0];
sx q[0];
rz(1.3657938) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.656865) q[2];
sx q[2];
rz(-0.51186168) q[2];
sx q[2];
rz(-0.32527015) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3615661) q[1];
sx q[1];
rz(-1.2793555) q[1];
sx q[1];
rz(2.7855278) q[1];
rz(-0.74541645) q[3];
sx q[3];
rz(-2.7099797) q[3];
sx q[3];
rz(3.0109757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9251755) q[2];
sx q[2];
rz(-1.5745796) q[2];
sx q[2];
rz(-0.65640059) q[2];
rz(2.9260855) q[3];
sx q[3];
rz(-1.7301205) q[3];
sx q[3];
rz(-1.6925252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0553174) q[0];
sx q[0];
rz(-1.3441939) q[0];
sx q[0];
rz(0.84981808) q[0];
rz(-0.96018106) q[1];
sx q[1];
rz(-2.7268867) q[1];
sx q[1];
rz(-0.93920952) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67341833) q[0];
sx q[0];
rz(-0.8020173) q[0];
sx q[0];
rz(1.9413663) q[0];
rz(-pi) q[1];
rz(1.1562111) q[2];
sx q[2];
rz(-1.5778981) q[2];
sx q[2];
rz(0.19609253) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85921849) q[1];
sx q[1];
rz(-1.9422429) q[1];
sx q[1];
rz(2.4747495) q[1];
rz(-pi) q[2];
rz(2.2978503) q[3];
sx q[3];
rz(-0.66202032) q[3];
sx q[3];
rz(-0.5211422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5822997) q[2];
sx q[2];
rz(-1.7377995) q[2];
sx q[2];
rz(-0.55029184) q[2];
rz(-3.0136717) q[3];
sx q[3];
rz(-3.0030799) q[3];
sx q[3];
rz(2.1730455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6561103) q[0];
sx q[0];
rz(-0.67018296) q[0];
sx q[0];
rz(-0.67076587) q[0];
rz(-2.7863964) q[1];
sx q[1];
rz(-1.6658446) q[1];
sx q[1];
rz(2.7459941) q[1];
rz(-2.998299) q[2];
sx q[2];
rz(-1.7096277) q[2];
sx q[2];
rz(-1.1451677) q[2];
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
