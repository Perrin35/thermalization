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
rz(0.056869153) q[0];
sx q[0];
rz(-0.19357227) q[0];
sx q[0];
rz(-0.40785664) q[0];
rz(-0.017539311) q[1];
sx q[1];
rz(-1.8899625) q[1];
sx q[1];
rz(-1.6012021) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9303794) q[0];
sx q[0];
rz(-1.5222852) q[0];
sx q[0];
rz(2.0082974) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45405252) q[2];
sx q[2];
rz(-1.6925998) q[2];
sx q[2];
rz(-1.6495088) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0527264) q[1];
sx q[1];
rz(-1.3590993) q[1];
sx q[1];
rz(-2.746341) q[1];
rz(-pi) q[2];
rz(-0.92272051) q[3];
sx q[3];
rz(-0.86258234) q[3];
sx q[3];
rz(-0.012383637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5876329) q[2];
sx q[2];
rz(-1.0593869) q[2];
sx q[2];
rz(-2.9239192) q[2];
rz(2.830128) q[3];
sx q[3];
rz(-2.5296827) q[3];
sx q[3];
rz(1.9757087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.72164732) q[0];
sx q[0];
rz(-2.8657275) q[0];
sx q[0];
rz(0.69197792) q[0];
rz(-0.75633374) q[1];
sx q[1];
rz(-2.6192009) q[1];
sx q[1];
rz(-1.5914894) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0333918) q[0];
sx q[0];
rz(-2.4188359) q[0];
sx q[0];
rz(0.18527822) q[0];
x q[1];
rz(1.9687551) q[2];
sx q[2];
rz(-1.6670456) q[2];
sx q[2];
rz(-0.41589662) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.34775388) q[1];
sx q[1];
rz(-1.4112541) q[1];
sx q[1];
rz(-1.6942339) q[1];
rz(-0.49868985) q[3];
sx q[3];
rz(-1.0674006) q[3];
sx q[3];
rz(1.4273912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1456566) q[2];
sx q[2];
rz(-0.17309509) q[2];
sx q[2];
rz(-0.23925979) q[2];
rz(-1.650882) q[3];
sx q[3];
rz(-1.2942634) q[3];
sx q[3];
rz(1.9132805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2759129) q[0];
sx q[0];
rz(-0.90621197) q[0];
sx q[0];
rz(-3.1269585) q[0];
rz(-0.89302653) q[1];
sx q[1];
rz(-0.97517401) q[1];
sx q[1];
rz(-0.21569529) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94616643) q[0];
sx q[0];
rz(-3.0901244) q[0];
sx q[0];
rz(1.6610751) q[0];
x q[1];
rz(1.9142308) q[2];
sx q[2];
rz(-0.94007713) q[2];
sx q[2];
rz(-1.3892) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1163997) q[1];
sx q[1];
rz(-1.4433675) q[1];
sx q[1];
rz(2.9078632) q[1];
rz(-2.2641597) q[3];
sx q[3];
rz(-1.0022396) q[3];
sx q[3];
rz(1.7753778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0821685) q[2];
sx q[2];
rz(-1.1801722) q[2];
sx q[2];
rz(2.1480985) q[2];
rz(0.70837402) q[3];
sx q[3];
rz(-2.0230484) q[3];
sx q[3];
rz(-0.12349252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74743903) q[0];
sx q[0];
rz(-0.49462947) q[0];
sx q[0];
rz(2.4476449) q[0];
rz(0.28678647) q[1];
sx q[1];
rz(-1.5319805) q[1];
sx q[1];
rz(-1.087711) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74831731) q[0];
sx q[0];
rz(-1.4469742) q[0];
sx q[0];
rz(-1.0626777) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55107815) q[2];
sx q[2];
rz(-1.113453) q[2];
sx q[2];
rz(-1.5975022) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3301311) q[1];
sx q[1];
rz(-2.4566659) q[1];
sx q[1];
rz(-1.0792988) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4272653) q[3];
sx q[3];
rz(-1.5363524) q[3];
sx q[3];
rz(-1.2726651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22471681) q[2];
sx q[2];
rz(-1.9852873) q[2];
sx q[2];
rz(1.3678331) q[2];
rz(-1.9035089) q[3];
sx q[3];
rz(-1.5725458) q[3];
sx q[3];
rz(0.21627538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5424159) q[0];
sx q[0];
rz(-1.2681862) q[0];
sx q[0];
rz(2.5237778) q[0];
rz(2.0333596) q[1];
sx q[1];
rz(-2.8384659) q[1];
sx q[1];
rz(2.2584426) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1681447) q[0];
sx q[0];
rz(-2.3855378) q[0];
sx q[0];
rz(0.010764695) q[0];
rz(-pi) q[1];
rz(-0.3381392) q[2];
sx q[2];
rz(-0.45881841) q[2];
sx q[2];
rz(0.79566075) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1146436) q[1];
sx q[1];
rz(-1.0806688) q[1];
sx q[1];
rz(0.074322083) q[1];
rz(-1.0155025) q[3];
sx q[3];
rz(-1.2153373) q[3];
sx q[3];
rz(-0.36830869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21657476) q[2];
sx q[2];
rz(-0.25187945) q[2];
sx q[2];
rz(1.3612755) q[2];
rz(-1.2733634) q[3];
sx q[3];
rz(-1.7073771) q[3];
sx q[3];
rz(-2.1586965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0971766) q[0];
sx q[0];
rz(-1.8955078) q[0];
sx q[0];
rz(-0.62028766) q[0];
rz(-2.7445131) q[1];
sx q[1];
rz(-2.670791) q[1];
sx q[1];
rz(1.9420067) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9387377) q[0];
sx q[0];
rz(-1.5498383) q[0];
sx q[0];
rz(1.6771132) q[0];
rz(-pi) q[1];
rz(-2.1753005) q[2];
sx q[2];
rz(-0.25560616) q[2];
sx q[2];
rz(2.5641455) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7677557) q[1];
sx q[1];
rz(-1.3779614) q[1];
sx q[1];
rz(2.6111433) q[1];
x q[2];
rz(-1.689246) q[3];
sx q[3];
rz(-1.9129686) q[3];
sx q[3];
rz(0.3698879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24641307) q[2];
sx q[2];
rz(-1.2440888) q[2];
sx q[2];
rz(0.67071521) q[2];
rz(2.8504168) q[3];
sx q[3];
rz(-2.5305735) q[3];
sx q[3];
rz(2.5197855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7156242) q[0];
sx q[0];
rz(-1.427303) q[0];
sx q[0];
rz(0.17844644) q[0];
rz(0.87002358) q[1];
sx q[1];
rz(-0.88302892) q[1];
sx q[1];
rz(-2.3387486) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74551302) q[0];
sx q[0];
rz(-2.1901696) q[0];
sx q[0];
rz(1.5399163) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58391352) q[2];
sx q[2];
rz(-1.7758992) q[2];
sx q[2];
rz(2.5896304) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4124914) q[1];
sx q[1];
rz(-1.8961398) q[1];
sx q[1];
rz(-2.4304683) q[1];
rz(-pi) q[2];
rz(2.1986897) q[3];
sx q[3];
rz(-1.3602837) q[3];
sx q[3];
rz(2.877416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37658438) q[2];
sx q[2];
rz(-1.8439801) q[2];
sx q[2];
rz(0.89361781) q[2];
rz(-1.5997959) q[3];
sx q[3];
rz(-1.4897646) q[3];
sx q[3];
rz(2.5347575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97745085) q[0];
sx q[0];
rz(-0.41005382) q[0];
sx q[0];
rz(-2.9264911) q[0];
rz(0.84456259) q[1];
sx q[1];
rz(-1.8212049) q[1];
sx q[1];
rz(-0.79427687) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0117019) q[0];
sx q[0];
rz(-0.029057682) q[0];
sx q[0];
rz(0.41883166) q[0];
rz(-0.5735917) q[2];
sx q[2];
rz(-2.6124138) q[2];
sx q[2];
rz(2.0647788) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8848443) q[1];
sx q[1];
rz(-1.0326385) q[1];
sx q[1];
rz(-3.0358008) q[1];
rz(-pi) q[2];
rz(-1.896554) q[3];
sx q[3];
rz(-0.81607258) q[3];
sx q[3];
rz(-1.2118424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8636785) q[2];
sx q[2];
rz(-1.9387551) q[2];
sx q[2];
rz(1.7318783) q[2];
rz(2.388741) q[3];
sx q[3];
rz(-1.9949621) q[3];
sx q[3];
rz(1.0021771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3442605) q[0];
sx q[0];
rz(-2.7334038) q[0];
sx q[0];
rz(2.1687188) q[0];
rz(-1.6196039) q[1];
sx q[1];
rz(-1.3643967) q[1];
sx q[1];
rz(0.63046986) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3909633) q[0];
sx q[0];
rz(-1.4328151) q[0];
sx q[0];
rz(0.18355455) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6182329) q[2];
sx q[2];
rz(-0.15711297) q[2];
sx q[2];
rz(-2.4542798) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4293824) q[1];
sx q[1];
rz(-2.3377067) q[1];
sx q[1];
rz(-1.6401119) q[1];
rz(-pi) q[2];
rz(1.6338324) q[3];
sx q[3];
rz(-1.3975289) q[3];
sx q[3];
rz(1.2931371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1260599) q[2];
sx q[2];
rz(-1.9893179) q[2];
sx q[2];
rz(-1.8401592) q[2];
rz(1.556501) q[3];
sx q[3];
rz(-2.3157178) q[3];
sx q[3];
rz(-2.7263156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62366098) q[0];
sx q[0];
rz(-2.9908337) q[0];
sx q[0];
rz(-2.6483722) q[0];
rz(-1.8863691) q[1];
sx q[1];
rz(-1.247765) q[1];
sx q[1];
rz(-2.7769322) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5887774) q[0];
sx q[0];
rz(-1.7087666) q[0];
sx q[0];
rz(1.2347883) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0027867) q[2];
sx q[2];
rz(-0.73633654) q[2];
sx q[2];
rz(-1.9970837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2233878) q[1];
sx q[1];
rz(-0.80092309) q[1];
sx q[1];
rz(0.078322874) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2544028) q[3];
sx q[3];
rz(-1.6091765) q[3];
sx q[3];
rz(0.61172133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6898592) q[2];
sx q[2];
rz(-1.006459) q[2];
sx q[2];
rz(-2.149392) q[2];
rz(2.774488) q[3];
sx q[3];
rz(-1.7016943) q[3];
sx q[3];
rz(2.4632857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51919666) q[0];
sx q[0];
rz(-1.665103) q[0];
sx q[0];
rz(-1.7892224) q[0];
rz(-2.2182111) q[1];
sx q[1];
rz(-0.8538178) q[1];
sx q[1];
rz(-1.1710844) q[1];
rz(-1.7708764) q[2];
sx q[2];
rz(-1.4309819) q[2];
sx q[2];
rz(-0.33071721) q[2];
rz(2.4685341) q[3];
sx q[3];
rz(-1.0792427) q[3];
sx q[3];
rz(-1.4918809) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
