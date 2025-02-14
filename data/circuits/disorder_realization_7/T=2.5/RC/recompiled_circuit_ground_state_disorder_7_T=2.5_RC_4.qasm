OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0291542) q[0];
sx q[0];
rz(-1.6310863) q[0];
sx q[0];
rz(-0.36614585) q[0];
rz(-1.0547628) q[1];
sx q[1];
rz(-1.2843479) q[1];
sx q[1];
rz(0.73959124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407362) q[0];
sx q[0];
rz(-1.2175517) q[0];
sx q[0];
rz(-2.059444) q[0];
x q[1];
rz(0.97402699) q[2];
sx q[2];
rz(-2.1197017) q[2];
sx q[2];
rz(3.10534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.057159) q[1];
sx q[1];
rz(-2.4600303) q[1];
sx q[1];
rz(0.8820028) q[1];
rz(-pi) q[2];
rz(-0.32259703) q[3];
sx q[3];
rz(-1.8290213) q[3];
sx q[3];
rz(0.74479529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71020484) q[2];
sx q[2];
rz(-1.7040161) q[2];
sx q[2];
rz(-2.0598742) q[2];
rz(-1.9234575) q[3];
sx q[3];
rz(-1.5798605) q[3];
sx q[3];
rz(-1.6703828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1350988) q[0];
sx q[0];
rz(-0.83469892) q[0];
sx q[0];
rz(2.304402) q[0];
rz(2.2039425) q[1];
sx q[1];
rz(-1.7559914) q[1];
sx q[1];
rz(-3.0167276) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0406348) q[0];
sx q[0];
rz(-2.0748108) q[0];
sx q[0];
rz(-2.0547778) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6972448) q[2];
sx q[2];
rz(-2.5160061) q[2];
sx q[2];
rz(1.7162965) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8229577) q[1];
sx q[1];
rz(-0.97170107) q[1];
sx q[1];
rz(1.4692719) q[1];
rz(-pi) q[2];
rz(0.30571823) q[3];
sx q[3];
rz(-1.5366712) q[3];
sx q[3];
rz(0.21473884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3426334) q[2];
sx q[2];
rz(-2.0755167) q[2];
sx q[2];
rz(2.7349045) q[2];
rz(2.0090571) q[3];
sx q[3];
rz(-1.5187902) q[3];
sx q[3];
rz(2.2741545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36658528) q[0];
sx q[0];
rz(-2.5560515) q[0];
sx q[0];
rz(-1.0803692) q[0];
rz(0.22251546) q[1];
sx q[1];
rz(-2.3718926) q[1];
sx q[1];
rz(2.5710107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2704318) q[0];
sx q[0];
rz(-0.50073114) q[0];
sx q[0];
rz(0.88520925) q[0];
rz(-pi) q[1];
rz(0.9502842) q[2];
sx q[2];
rz(-2.0740905) q[2];
sx q[2];
rz(2.2719943) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1751082) q[1];
sx q[1];
rz(-2.3382332) q[1];
sx q[1];
rz(0.030637189) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57502745) q[3];
sx q[3];
rz(-1.417727) q[3];
sx q[3];
rz(2.0401772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.83501619) q[2];
sx q[2];
rz(-1.4811652) q[2];
sx q[2];
rz(1.906685) q[2];
rz(-2.0645449) q[3];
sx q[3];
rz(-1.5826694) q[3];
sx q[3];
rz(-0.41729116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1398337) q[0];
sx q[0];
rz(-1.3707021) q[0];
sx q[0];
rz(2.4613001) q[0];
rz(1.3975877) q[1];
sx q[1];
rz(-2.0060507) q[1];
sx q[1];
rz(0.23522338) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096073791) q[0];
sx q[0];
rz(-2.1102724) q[0];
sx q[0];
rz(1.1474613) q[0];
x q[1];
rz(0.3346457) q[2];
sx q[2];
rz(-0.88289936) q[2];
sx q[2];
rz(-2.2918224) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1295812) q[1];
sx q[1];
rz(-2.5223603) q[1];
sx q[1];
rz(2.7577452) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8743478) q[3];
sx q[3];
rz(-1.6755723) q[3];
sx q[3];
rz(-1.1997833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.32517165) q[2];
sx q[2];
rz(-2.0656526) q[2];
sx q[2];
rz(-1.761033) q[2];
rz(-2.374968) q[3];
sx q[3];
rz(-2.4216757) q[3];
sx q[3];
rz(0.64002526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0470873) q[0];
sx q[0];
rz(-2.6386059) q[0];
sx q[0];
rz(-1.267953) q[0];
rz(-0.9371593) q[1];
sx q[1];
rz(-2.2068534) q[1];
sx q[1];
rz(2.8840205) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99058305) q[0];
sx q[0];
rz(-2.5611612) q[0];
sx q[0];
rz(1.2925757) q[0];
rz(2.0054162) q[2];
sx q[2];
rz(-2.4939257) q[2];
sx q[2];
rz(2.2369564) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7422545) q[1];
sx q[1];
rz(-0.4961001) q[1];
sx q[1];
rz(0.45544405) q[1];
rz(-pi) q[2];
rz(-2.0553642) q[3];
sx q[3];
rz(-1.7811097) q[3];
sx q[3];
rz(-2.9919101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3246954) q[2];
sx q[2];
rz(-2.2163053) q[2];
sx q[2];
rz(2.2980105) q[2];
rz(2.2364565) q[3];
sx q[3];
rz(-2.2456808) q[3];
sx q[3];
rz(2.8781273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0819241) q[0];
sx q[0];
rz(-0.34316871) q[0];
sx q[0];
rz(-1.7053509) q[0];
rz(-1.2417271) q[1];
sx q[1];
rz(-1.4271586) q[1];
sx q[1];
rz(2.5232975) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70253187) q[0];
sx q[0];
rz(-1.0510604) q[0];
sx q[0];
rz(2.5072839) q[0];
rz(-pi) q[1];
rz(-1.3682057) q[2];
sx q[2];
rz(-2.4491022) q[2];
sx q[2];
rz(-0.67815526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8346295) q[1];
sx q[1];
rz(-1.2696046) q[1];
sx q[1];
rz(-2.7121915) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6446436) q[3];
sx q[3];
rz(-2.5922311) q[3];
sx q[3];
rz(-0.13298154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5228086) q[2];
sx q[2];
rz(-1.2455384) q[2];
sx q[2];
rz(-1.9469384) q[2];
rz(1.4905802) q[3];
sx q[3];
rz(-1.4158764) q[3];
sx q[3];
rz(1.850261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2456583) q[0];
sx q[0];
rz(-2.7722562) q[0];
sx q[0];
rz(-2.9781407) q[0];
rz(-2.5864736) q[1];
sx q[1];
rz(-1.6682383) q[1];
sx q[1];
rz(-0.82383531) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6763944) q[0];
sx q[0];
rz(-1.4270743) q[0];
sx q[0];
rz(0.31938796) q[0];
rz(-0.78764913) q[2];
sx q[2];
rz(-1.3795492) q[2];
sx q[2];
rz(2.8764962) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2545611) q[1];
sx q[1];
rz(-1.3929092) q[1];
sx q[1];
rz(-0.46079854) q[1];
rz(-2.7300445) q[3];
sx q[3];
rz(-1.8754861) q[3];
sx q[3];
rz(2.2652486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0550363) q[2];
sx q[2];
rz(-1.8359416) q[2];
sx q[2];
rz(1.6769064) q[2];
rz(-0.054051789) q[3];
sx q[3];
rz(-0.87658221) q[3];
sx q[3];
rz(-2.7189253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61580491) q[0];
sx q[0];
rz(-2.7532888) q[0];
sx q[0];
rz(1.4317321) q[0];
rz(1.8852385) q[1];
sx q[1];
rz(-1.6836555) q[1];
sx q[1];
rz(1.8690522) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1521069) q[0];
sx q[0];
rz(-2.4740334) q[0];
sx q[0];
rz(2.9015673) q[0];
rz(1.854004) q[2];
sx q[2];
rz(-2.0618366) q[2];
sx q[2];
rz(-1.1591737) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1450729) q[1];
sx q[1];
rz(-0.47059083) q[1];
sx q[1];
rz(0.7820635) q[1];
x q[2];
rz(-0.36964396) q[3];
sx q[3];
rz(-2.3497006) q[3];
sx q[3];
rz(1.7580851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.58497477) q[2];
sx q[2];
rz(-1.40404) q[2];
sx q[2];
rz(-0.72648826) q[2];
rz(0.11897421) q[3];
sx q[3];
rz(-1.3930895) q[3];
sx q[3];
rz(0.52072853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.2105836) q[0];
sx q[0];
rz(-1.9887661) q[0];
sx q[0];
rz(-0.3535122) q[0];
rz(1.1734236) q[1];
sx q[1];
rz(-1.8559772) q[1];
sx q[1];
rz(-1.2849503) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.374732) q[0];
sx q[0];
rz(-2.6832125) q[0];
sx q[0];
rz(-1.2552257) q[0];
rz(-2.6899723) q[2];
sx q[2];
rz(-1.2112144) q[2];
sx q[2];
rz(3.0082085) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.28390005) q[1];
sx q[1];
rz(-1.5491747) q[1];
sx q[1];
rz(1.2766854) q[1];
rz(-pi) q[2];
rz(-2.5250149) q[3];
sx q[3];
rz(-1.7093466) q[3];
sx q[3];
rz(-1.5091173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.53359199) q[2];
sx q[2];
rz(-2.3054391) q[2];
sx q[2];
rz(2.3261435) q[2];
rz(2.4943374) q[3];
sx q[3];
rz(-1.7630968) q[3];
sx q[3];
rz(0.7816202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4746998) q[0];
sx q[0];
rz(-0.26645461) q[0];
sx q[0];
rz(1.4315963) q[0];
rz(-2.73009) q[1];
sx q[1];
rz(-1.9930379) q[1];
sx q[1];
rz(3.0130951) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.381739) q[0];
sx q[0];
rz(-2.2138192) q[0];
sx q[0];
rz(-1.169315) q[0];
rz(3.0017486) q[2];
sx q[2];
rz(-2.2763414) q[2];
sx q[2];
rz(-1.2024513) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.45396475) q[1];
sx q[1];
rz(-1.0610136) q[1];
sx q[1];
rz(2.3483454) q[1];
rz(0.64537489) q[3];
sx q[3];
rz(-2.5119492) q[3];
sx q[3];
rz(-0.30984344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0633885) q[2];
sx q[2];
rz(-1.3271164) q[2];
sx q[2];
rz(-0.77888954) q[2];
rz(-1.8552467) q[3];
sx q[3];
rz(-2.0467919) q[3];
sx q[3];
rz(-1.8761427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.417199) q[0];
sx q[0];
rz(-1.9211641) q[0];
sx q[0];
rz(-1.7838508) q[0];
rz(2.4667274) q[1];
sx q[1];
rz(-1.6541031) q[1];
sx q[1];
rz(-2.0872603) q[1];
rz(0.40140006) q[2];
sx q[2];
rz(-1.4742322) q[2];
sx q[2];
rz(-0.36591224) q[2];
rz(-1.5065083) q[3];
sx q[3];
rz(-2.3953606) q[3];
sx q[3];
rz(2.2904979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
