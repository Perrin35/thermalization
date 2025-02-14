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
rz(-0.80225575) q[0];
sx q[0];
rz(-1.7576317) q[0];
sx q[0];
rz(1.2686165) q[0];
rz(0.4624548) q[1];
sx q[1];
rz(-3.0825244) q[1];
sx q[1];
rz(1.4055835) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2911513) q[0];
sx q[0];
rz(-1.1014525) q[0];
sx q[0];
rz(-0.48908451) q[0];
rz(-0.67063801) q[2];
sx q[2];
rz(-1.3893638) q[2];
sx q[2];
rz(1.8701613) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51912145) q[1];
sx q[1];
rz(-2.2833385) q[1];
sx q[1];
rz(2.126241) q[1];
rz(-2.2470993) q[3];
sx q[3];
rz(-1.749012) q[3];
sx q[3];
rz(-1.3823929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1543701) q[2];
sx q[2];
rz(-1.9742249) q[2];
sx q[2];
rz(-0.78987375) q[2];
rz(3.1176944) q[3];
sx q[3];
rz(-1.8022715) q[3];
sx q[3];
rz(-2.6051615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24870482) q[0];
sx q[0];
rz(-1.4905812) q[0];
sx q[0];
rz(1.254068) q[0];
rz(1.4887384) q[1];
sx q[1];
rz(-2.2658927) q[1];
sx q[1];
rz(-0.48315963) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7308759) q[0];
sx q[0];
rz(-2.9456101) q[0];
sx q[0];
rz(1.3021126) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4065521) q[2];
sx q[2];
rz(-1.9695821) q[2];
sx q[2];
rz(0.48828188) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4401958) q[1];
sx q[1];
rz(-0.45397511) q[1];
sx q[1];
rz(-2.5099436) q[1];
x q[2];
rz(-1.3951357) q[3];
sx q[3];
rz(-1.2845728) q[3];
sx q[3];
rz(1.7475048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2529605) q[2];
sx q[2];
rz(-1.4560207) q[2];
sx q[2];
rz(-1.6271094) q[2];
rz(3.0857981) q[3];
sx q[3];
rz(-1.5972219) q[3];
sx q[3];
rz(1.4896711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91419879) q[0];
sx q[0];
rz(-2.6970503) q[0];
sx q[0];
rz(0.12810531) q[0];
rz(-0.57614342) q[1];
sx q[1];
rz(-2.4100401) q[1];
sx q[1];
rz(2.3760956) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.573581) q[0];
sx q[0];
rz(-0.67493248) q[0];
sx q[0];
rz(0.57008596) q[0];
rz(-pi) q[1];
rz(0.68674318) q[2];
sx q[2];
rz(-1.8675065) q[2];
sx q[2];
rz(-2.3627757) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0772503) q[1];
sx q[1];
rz(-2.0301452) q[1];
sx q[1];
rz(0.5456173) q[1];
x q[2];
rz(-2.9981444) q[3];
sx q[3];
rz(-1.5455523) q[3];
sx q[3];
rz(-0.48893602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8219882) q[2];
sx q[2];
rz(-0.87279785) q[2];
sx q[2];
rz(-2.2785211) q[2];
rz(-2.7989164) q[3];
sx q[3];
rz(-0.42049146) q[3];
sx q[3];
rz(-1.4884523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89219013) q[0];
sx q[0];
rz(-0.96977314) q[0];
sx q[0];
rz(-0.48200193) q[0];
rz(-2.8328698) q[1];
sx q[1];
rz(-0.32544193) q[1];
sx q[1];
rz(-0.6122922) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68744874) q[0];
sx q[0];
rz(-0.73546919) q[0];
sx q[0];
rz(-0.81540458) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6545534) q[2];
sx q[2];
rz(-1.5568019) q[2];
sx q[2];
rz(2.9521041) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6867076) q[1];
sx q[1];
rz(-2.5338182) q[1];
sx q[1];
rz(-1.6493504) q[1];
rz(-pi) q[2];
rz(0.14928603) q[3];
sx q[3];
rz(-1.5773492) q[3];
sx q[3];
rz(-0.1803151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4838532) q[2];
sx q[2];
rz(-2.190399) q[2];
sx q[2];
rz(-1.308002) q[2];
rz(0.43404239) q[3];
sx q[3];
rz(-1.3515892) q[3];
sx q[3];
rz(1.8724117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16172116) q[0];
sx q[0];
rz(-1.5890108) q[0];
sx q[0];
rz(2.8671434) q[0];
rz(1.5212003) q[1];
sx q[1];
rz(-2.0976286) q[1];
sx q[1];
rz(1.257198) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6406157) q[0];
sx q[0];
rz(-1.8145992) q[0];
sx q[0];
rz(-1.8200726) q[0];
x q[1];
rz(2.1047701) q[2];
sx q[2];
rz(-1.8920915) q[2];
sx q[2];
rz(0.23201135) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8179743) q[1];
sx q[1];
rz(-2.5516641) q[1];
sx q[1];
rz(1.1765615) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86602314) q[3];
sx q[3];
rz(-1.3228205) q[3];
sx q[3];
rz(-3.0803404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3955128) q[2];
sx q[2];
rz(-1.6359676) q[2];
sx q[2];
rz(-2.5858322) q[2];
rz(2.3505576) q[3];
sx q[3];
rz(-1.0020703) q[3];
sx q[3];
rz(1.6382431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24285862) q[0];
sx q[0];
rz(-1.7019615) q[0];
sx q[0];
rz(-2.4295501) q[0];
rz(2.5391319) q[1];
sx q[1];
rz(-0.41350499) q[1];
sx q[1];
rz(0.079040225) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1940546) q[0];
sx q[0];
rz(-1.6350766) q[0];
sx q[0];
rz(0.044919515) q[0];
rz(-pi) q[1];
rz(-1.4858133) q[2];
sx q[2];
rz(-1.6459609) q[2];
sx q[2];
rz(-2.4369881) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1823719) q[1];
sx q[1];
rz(-0.8405662) q[1];
sx q[1];
rz(-1.8664136) q[1];
x q[2];
rz(1.4677804) q[3];
sx q[3];
rz(-1.2742189) q[3];
sx q[3];
rz(-1.0271343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1003803) q[2];
sx q[2];
rz(-2.2104287) q[2];
sx q[2];
rz(2.1616914) q[2];
rz(1.6117217) q[3];
sx q[3];
rz(-1.2414705) q[3];
sx q[3];
rz(2.8821168) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74705446) q[0];
sx q[0];
rz(-1.2815481) q[0];
sx q[0];
rz(3.0406612) q[0];
rz(-2.586567) q[1];
sx q[1];
rz(-0.9340159) q[1];
sx q[1];
rz(1.2947882) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053558401) q[0];
sx q[0];
rz(-0.71672601) q[0];
sx q[0];
rz(1.929856) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1406222) q[2];
sx q[2];
rz(-1.8651267) q[2];
sx q[2];
rz(-0.4145588) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6375512) q[1];
sx q[1];
rz(-2.6890272) q[1];
sx q[1];
rz(-0.33337793) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0247308) q[3];
sx q[3];
rz(-2.6203558) q[3];
sx q[3];
rz(1.8841528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1780221) q[2];
sx q[2];
rz(-2.2691085) q[2];
sx q[2];
rz(-2.8998609) q[2];
rz(2.669615) q[3];
sx q[3];
rz(-0.21502544) q[3];
sx q[3];
rz(1.5579461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3374775) q[0];
sx q[0];
rz(-0.98144704) q[0];
sx q[0];
rz(-0.61019439) q[0];
rz(-0.21774165) q[1];
sx q[1];
rz(-0.95548958) q[1];
sx q[1];
rz(-2.311923) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.037704) q[0];
sx q[0];
rz(-2.12426) q[0];
sx q[0];
rz(0.82413701) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83135624) q[2];
sx q[2];
rz(-2.6717253) q[2];
sx q[2];
rz(2.510315) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9431994) q[1];
sx q[1];
rz(-1.7560609) q[1];
sx q[1];
rz(-2.7279502) q[1];
rz(2.2822718) q[3];
sx q[3];
rz(-1.4830657) q[3];
sx q[3];
rz(2.6130207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8391476) q[2];
sx q[2];
rz(-1.382901) q[2];
sx q[2];
rz(1.1350606) q[2];
rz(-2.6531687) q[3];
sx q[3];
rz(-1.6596183) q[3];
sx q[3];
rz(-1.9058913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9957073) q[0];
sx q[0];
rz(-1.1111525) q[0];
sx q[0];
rz(-2.1671894) q[0];
rz(0.79565945) q[1];
sx q[1];
rz(-2.7641422) q[1];
sx q[1];
rz(0.61765751) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5461676) q[0];
sx q[0];
rz(-1.6542395) q[0];
sx q[0];
rz(1.5682778) q[0];
x q[1];
rz(1.9699924) q[2];
sx q[2];
rz(-2.0621057) q[2];
sx q[2];
rz(2.0727061) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61231919) q[1];
sx q[1];
rz(-2.5041083) q[1];
sx q[1];
rz(-0.46302621) q[1];
rz(0.32989721) q[3];
sx q[3];
rz(-1.1217692) q[3];
sx q[3];
rz(-0.87775081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.016102942) q[2];
sx q[2];
rz(-1.8718655) q[2];
sx q[2];
rz(-1.3886836) q[2];
rz(-1.3056508) q[3];
sx q[3];
rz(-0.7235705) q[3];
sx q[3];
rz(-1.8587662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1284803) q[0];
sx q[0];
rz(-0.73902577) q[0];
sx q[0];
rz(-0.59984961) q[0];
rz(-2.5685617) q[1];
sx q[1];
rz(-0.60879469) q[1];
sx q[1];
rz(-2.1175687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1896418) q[0];
sx q[0];
rz(-1.1770403) q[0];
sx q[0];
rz(-2.0102565) q[0];
x q[1];
rz(1.5285792) q[2];
sx q[2];
rz(-0.91860549) q[2];
sx q[2];
rz(0.025328115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9309195) q[1];
sx q[1];
rz(-1.1751047) q[1];
sx q[1];
rz(1.4798505) q[1];
x q[2];
rz(1.7217595) q[3];
sx q[3];
rz(-1.1929885) q[3];
sx q[3];
rz(-2.092416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92274252) q[2];
sx q[2];
rz(-1.3877733) q[2];
sx q[2];
rz(2.3663523) q[2];
rz(1.6058589) q[3];
sx q[3];
rz(-1.1865059) q[3];
sx q[3];
rz(-1.8839914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(2.6105462) q[0];
sx q[0];
rz(-2.1903867) q[0];
sx q[0];
rz(1.2280986) q[0];
rz(-1.7465406) q[1];
sx q[1];
rz(-1.6680622) q[1];
sx q[1];
rz(-0.39573085) q[1];
rz(-1.878123) q[2];
sx q[2];
rz(-1.1358481) q[2];
sx q[2];
rz(-0.67107558) q[2];
rz(-0.80202924) q[3];
sx q[3];
rz(-1.5939972) q[3];
sx q[3];
rz(0.93929285) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
