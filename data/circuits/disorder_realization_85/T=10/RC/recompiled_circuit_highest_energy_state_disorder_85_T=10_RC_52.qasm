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
rz(-2.149481) q[0];
sx q[0];
rz(-0.7465201) q[0];
sx q[0];
rz(2.574805) q[0];
rz(1.2504638) q[1];
sx q[1];
rz(-1.992978) q[1];
sx q[1];
rz(-0.9196547) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84595478) q[0];
sx q[0];
rz(-1.7311117) q[0];
sx q[0];
rz(-0.79248595) q[0];
rz(-pi) q[1];
rz(2.5819725) q[2];
sx q[2];
rz(-2.3336448) q[2];
sx q[2];
rz(-1.9576278) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30324591) q[1];
sx q[1];
rz(-0.64323264) q[1];
sx q[1];
rz(-3.0641765) q[1];
rz(-pi) q[2];
rz(-2.9018974) q[3];
sx q[3];
rz(-2.5928232) q[3];
sx q[3];
rz(-1.5645998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2016242) q[2];
sx q[2];
rz(-2.2087966) q[2];
sx q[2];
rz(1.0533818) q[2];
rz(0.71422226) q[3];
sx q[3];
rz(-0.37319365) q[3];
sx q[3];
rz(-1.8478954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6222222) q[0];
sx q[0];
rz(-0.70835963) q[0];
sx q[0];
rz(-0.57193065) q[0];
rz(-1.8427303) q[1];
sx q[1];
rz(-0.79890257) q[1];
sx q[1];
rz(0.78972185) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.84452) q[0];
sx q[0];
rz(-1.6668586) q[0];
sx q[0];
rz(1.8734832) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65424325) q[2];
sx q[2];
rz(-1.6983319) q[2];
sx q[2];
rz(-0.5687643) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4724628) q[1];
sx q[1];
rz(-1.2910559) q[1];
sx q[1];
rz(2.521924) q[1];
rz(-pi) q[2];
rz(2.5511207) q[3];
sx q[3];
rz(-2.3749224) q[3];
sx q[3];
rz(1.9719415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.85084891) q[2];
sx q[2];
rz(-0.76883832) q[2];
sx q[2];
rz(1.3624462) q[2];
rz(0.29081523) q[3];
sx q[3];
rz(-1.3664061) q[3];
sx q[3];
rz(3.0349558) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0027851) q[0];
sx q[0];
rz(-1.0951575) q[0];
sx q[0];
rz(-0.7005257) q[0];
rz(-0.48577148) q[1];
sx q[1];
rz(-2.3120717) q[1];
sx q[1];
rz(-0.22451678) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.016708) q[0];
sx q[0];
rz(-1.9210235) q[0];
sx q[0];
rz(1.4136397) q[0];
rz(2.1491884) q[2];
sx q[2];
rz(-1.4881721) q[2];
sx q[2];
rz(-2.3985661) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2292994) q[1];
sx q[1];
rz(-2.5535431) q[1];
sx q[1];
rz(-1.4347784) q[1];
rz(-pi) q[2];
x q[2];
rz(1.443601) q[3];
sx q[3];
rz(-0.7326829) q[3];
sx q[3];
rz(-0.81058433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.756788) q[2];
sx q[2];
rz(-2.2497358) q[2];
sx q[2];
rz(-3.0925114) q[2];
rz(-3.0745506) q[3];
sx q[3];
rz(-1.9837572) q[3];
sx q[3];
rz(0.66655794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9662358) q[0];
sx q[0];
rz(-0.55713621) q[0];
sx q[0];
rz(-3.1224342) q[0];
rz(1.3592023) q[1];
sx q[1];
rz(-1.4298341) q[1];
sx q[1];
rz(-0.0044936831) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.27994) q[0];
sx q[0];
rz(-0.77095448) q[0];
sx q[0];
rz(2.3907135) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49764823) q[2];
sx q[2];
rz(-1.9714173) q[2];
sx q[2];
rz(1.7965574) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7178065) q[1];
sx q[1];
rz(-2.7423334) q[1];
sx q[1];
rz(-2.0171793) q[1];
x q[2];
rz(2.717797) q[3];
sx q[3];
rz(-2.9188699) q[3];
sx q[3];
rz(0.84171695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6898474) q[2];
sx q[2];
rz(-1.0901901) q[2];
sx q[2];
rz(-1.6881662) q[2];
rz(1.2545741) q[3];
sx q[3];
rz(-1.5213608) q[3];
sx q[3];
rz(0.5622676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8101863) q[0];
sx q[0];
rz(-1.9047381) q[0];
sx q[0];
rz(1.1619262) q[0];
rz(-2.3792073) q[1];
sx q[1];
rz(-0.98680174) q[1];
sx q[1];
rz(-0.42957482) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21107656) q[0];
sx q[0];
rz(-1.1628224) q[0];
sx q[0];
rz(0.039647722) q[0];
rz(-pi) q[1];
rz(1.8354906) q[2];
sx q[2];
rz(-2.3117723) q[2];
sx q[2];
rz(3.004937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1112238) q[1];
sx q[1];
rz(-2.3850523) q[1];
sx q[1];
rz(-1.0683505) q[1];
x q[2];
rz(2.0031092) q[3];
sx q[3];
rz(-0.70037445) q[3];
sx q[3];
rz(-2.5770503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9261711) q[2];
sx q[2];
rz(-1.2913707) q[2];
sx q[2];
rz(-1.8298836) q[2];
rz(2.6808776) q[3];
sx q[3];
rz(-2.1381133) q[3];
sx q[3];
rz(3.0084012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8120414) q[0];
sx q[0];
rz(-0.1596182) q[0];
sx q[0];
rz(1.106369) q[0];
rz(0.1144935) q[1];
sx q[1];
rz(-2.0888927) q[1];
sx q[1];
rz(0.053650275) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44941266) q[0];
sx q[0];
rz(-2.8859038) q[0];
sx q[0];
rz(1.5935672) q[0];
rz(-pi) q[1];
rz(0.099951115) q[2];
sx q[2];
rz(-2.3315773) q[2];
sx q[2];
rz(2.9423327) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70354706) q[1];
sx q[1];
rz(-1.6247735) q[1];
sx q[1];
rz(1.0590068) q[1];
x q[2];
rz(-2.5074282) q[3];
sx q[3];
rz(-1.0552707) q[3];
sx q[3];
rz(-0.17532119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8657118) q[2];
sx q[2];
rz(-1.2578473) q[2];
sx q[2];
rz(-2.6344521) q[2];
rz(1.7670613) q[3];
sx q[3];
rz(-1.5406698) q[3];
sx q[3];
rz(1.9929632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.62992612) q[0];
sx q[0];
rz(-2.0138795) q[0];
sx q[0];
rz(3.1003057) q[0];
rz(-0.73829007) q[1];
sx q[1];
rz(-1.2246882) q[1];
sx q[1];
rz(-0.98947492) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6302297) q[0];
sx q[0];
rz(-1.0697027) q[0];
sx q[0];
rz(2.5998235) q[0];
rz(-pi) q[1];
rz(-1.3812387) q[2];
sx q[2];
rz(-1.9206534) q[2];
sx q[2];
rz(2.7477392) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.38997) q[1];
sx q[1];
rz(-2.1643359) q[1];
sx q[1];
rz(-0.71243993) q[1];
x q[2];
rz(-1.0999849) q[3];
sx q[3];
rz(-0.71094497) q[3];
sx q[3];
rz(-0.31490745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5993293) q[2];
sx q[2];
rz(-1.0487391) q[2];
sx q[2];
rz(-0.85912022) q[2];
rz(-3.1298992) q[3];
sx q[3];
rz(-1.6418567) q[3];
sx q[3];
rz(-0.11277994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82666731) q[0];
sx q[0];
rz(-2.5496917) q[0];
sx q[0];
rz(0.23183204) q[0];
rz(-2.5595317) q[1];
sx q[1];
rz(-1.3287909) q[1];
sx q[1];
rz(1.9225165) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1506596) q[0];
sx q[0];
rz(-2.5703589) q[0];
sx q[0];
rz(-0.21749638) q[0];
rz(-pi) q[1];
rz(-3.0000971) q[2];
sx q[2];
rz(-2.7853051) q[2];
sx q[2];
rz(-2.3959999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.92182175) q[1];
sx q[1];
rz(-1.3857462) q[1];
sx q[1];
rz(-2.9875523) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0295415) q[3];
sx q[3];
rz(-0.89139056) q[3];
sx q[3];
rz(2.816659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9570738) q[2];
sx q[2];
rz(-1.3513214) q[2];
sx q[2];
rz(-0.78682023) q[2];
rz(1.5015073) q[3];
sx q[3];
rz(-0.4314751) q[3];
sx q[3];
rz(-1.4260346) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96104923) q[0];
sx q[0];
rz(-1.7683872) q[0];
sx q[0];
rz(-0.94026646) q[0];
rz(-0.44479784) q[1];
sx q[1];
rz(-2.2856789) q[1];
sx q[1];
rz(-0.14642265) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43035591) q[0];
sx q[0];
rz(-0.32546639) q[0];
sx q[0];
rz(-2.7379166) q[0];
rz(-pi) q[1];
rz(1.3166041) q[2];
sx q[2];
rz(-1.3066402) q[2];
sx q[2];
rz(3.11655) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6633353) q[1];
sx q[1];
rz(-0.90090776) q[1];
sx q[1];
rz(1.2557593) q[1];
rz(-3.0412263) q[3];
sx q[3];
rz(-2.3477738) q[3];
sx q[3];
rz(0.19089261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0545097) q[2];
sx q[2];
rz(-1.4887709) q[2];
sx q[2];
rz(-0.83756891) q[2];
rz(-2.2086823) q[3];
sx q[3];
rz(-2.1541903) q[3];
sx q[3];
rz(-2.3605409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8785716) q[0];
sx q[0];
rz(-1.913338) q[0];
sx q[0];
rz(-1.3680869) q[0];
rz(-0.076016501) q[1];
sx q[1];
rz(-0.58034211) q[1];
sx q[1];
rz(2.0215633) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0033542643) q[0];
sx q[0];
rz(-1.8347667) q[0];
sx q[0];
rz(1.2498943) q[0];
x q[1];
rz(-2.7897808) q[2];
sx q[2];
rz(-2.4016671) q[2];
sx q[2];
rz(-0.49525317) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8260086) q[1];
sx q[1];
rz(-2.8346229) q[1];
sx q[1];
rz(-1.2493285) q[1];
x q[2];
rz(1.4073652) q[3];
sx q[3];
rz(-1.4686958) q[3];
sx q[3];
rz(-0.40615505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40176216) q[2];
sx q[2];
rz(-1.6900475) q[2];
sx q[2];
rz(-0.24701992) q[2];
rz(0.57010993) q[3];
sx q[3];
rz(-2.6542122) q[3];
sx q[3];
rz(2.0232239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0433255) q[0];
sx q[0];
rz(-1.7541616) q[0];
sx q[0];
rz(2.7761205) q[0];
rz(3.0430766) q[1];
sx q[1];
rz(-2.5010074) q[1];
sx q[1];
rz(0.88190257) q[1];
rz(-1.1578887) q[2];
sx q[2];
rz(-1.5831686) q[2];
sx q[2];
rz(-1.8148212) q[2];
rz(2.5539342) q[3];
sx q[3];
rz(-0.44787221) q[3];
sx q[3];
rz(-0.032240562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
