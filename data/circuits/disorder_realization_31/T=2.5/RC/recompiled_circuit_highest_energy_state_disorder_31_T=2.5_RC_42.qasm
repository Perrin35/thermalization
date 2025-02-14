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
rz(0.54654044) q[0];
sx q[0];
rz(3.300907) q[0];
sx q[0];
rz(11.309639) q[0];
rz(0.15428421) q[1];
sx q[1];
rz(-2.0576539) q[1];
sx q[1];
rz(1.3561603) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.529146) q[0];
sx q[0];
rz(-0.52781289) q[0];
sx q[0];
rz(1.1124658) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9262395) q[2];
sx q[2];
rz(-0.15654187) q[2];
sx q[2];
rz(2.6555344) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.034223864) q[1];
sx q[1];
rz(-1.6315206) q[1];
sx q[1];
rz(0.71141457) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47812652) q[3];
sx q[3];
rz(-1.193913) q[3];
sx q[3];
rz(1.865975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.95947444) q[2];
sx q[2];
rz(-2.0017767) q[2];
sx q[2];
rz(-0.095890447) q[2];
rz(0.31283665) q[3];
sx q[3];
rz(-1.0597119) q[3];
sx q[3];
rz(-0.50725168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4023912) q[0];
sx q[0];
rz(-1.6283789) q[0];
sx q[0];
rz(-1.8550523) q[0];
rz(1.8402428) q[1];
sx q[1];
rz(-1.2570612) q[1];
sx q[1];
rz(1.4305065) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.032784) q[0];
sx q[0];
rz(-0.51575249) q[0];
sx q[0];
rz(-1.4620251) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3545121) q[2];
sx q[2];
rz(-2.5445017) q[2];
sx q[2];
rz(2.7379089) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7315361) q[1];
sx q[1];
rz(-1.441873) q[1];
sx q[1];
rz(-1.6881264) q[1];
x q[2];
rz(1.3271403) q[3];
sx q[3];
rz(-2.2639963) q[3];
sx q[3];
rz(1.5290631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8282738) q[2];
sx q[2];
rz(-2.2602849) q[2];
sx q[2];
rz(0.71448294) q[2];
rz(1.03164) q[3];
sx q[3];
rz(-2.1561626) q[3];
sx q[3];
rz(0.45903444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3752876) q[0];
sx q[0];
rz(-2.3630688) q[0];
sx q[0];
rz(0.2226204) q[0];
rz(-2.8011232) q[1];
sx q[1];
rz(-1.4687931) q[1];
sx q[1];
rz(0.31390831) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71843159) q[0];
sx q[0];
rz(-0.88215798) q[0];
sx q[0];
rz(-0.60899831) q[0];
rz(-pi) q[1];
rz(-1.0386499) q[2];
sx q[2];
rz(-2.3380816) q[2];
sx q[2];
rz(-1.6626128) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7463689) q[1];
sx q[1];
rz(-0.88942553) q[1];
sx q[1];
rz(-1.9699775) q[1];
rz(-pi) q[2];
rz(1.7944853) q[3];
sx q[3];
rz(-1.0850226) q[3];
sx q[3];
rz(0.94625326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8935304) q[2];
sx q[2];
rz(-1.6174822) q[2];
sx q[2];
rz(-2.1171872) q[2];
rz(-1.3192568) q[3];
sx q[3];
rz(-2.8596467) q[3];
sx q[3];
rz(-2.8446741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.1612741) q[0];
sx q[0];
rz(-1.7151105) q[0];
sx q[0];
rz(2.1153765) q[0];
rz(-0.16464344) q[1];
sx q[1];
rz(-2.7963729) q[1];
sx q[1];
rz(-0.77273291) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5467632) q[0];
sx q[0];
rz(-1.208361) q[0];
sx q[0];
rz(-0.67154638) q[0];
x q[1];
rz(1.5667138) q[2];
sx q[2];
rz(-1.2130623) q[2];
sx q[2];
rz(-0.28752354) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31269909) q[1];
sx q[1];
rz(-0.71694198) q[1];
sx q[1];
rz(-1.2374876) q[1];
rz(-pi) q[2];
rz(1.4882632) q[3];
sx q[3];
rz(-0.051710796) q[3];
sx q[3];
rz(-2.5178096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2356977) q[2];
sx q[2];
rz(-0.57196456) q[2];
sx q[2];
rz(-0.31615654) q[2];
rz(0.21137992) q[3];
sx q[3];
rz(-1.4485161) q[3];
sx q[3];
rz(-2.0785418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0036156) q[0];
sx q[0];
rz(-1.6162385) q[0];
sx q[0];
rz(2.3620918) q[0];
rz(-1.4315073) q[1];
sx q[1];
rz(-1.5305488) q[1];
sx q[1];
rz(3.0140108) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6955285) q[0];
sx q[0];
rz(-1.4565153) q[0];
sx q[0];
rz(-3.1067239) q[0];
rz(-pi) q[1];
x q[1];
rz(2.374619) q[2];
sx q[2];
rz(-0.80000814) q[2];
sx q[2];
rz(0.50398477) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0810839) q[1];
sx q[1];
rz(-1.2049872) q[1];
sx q[1];
rz(-2.8660956) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5041994) q[3];
sx q[3];
rz(-1.069832) q[3];
sx q[3];
rz(2.4911309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7369507) q[2];
sx q[2];
rz(-1.5664132) q[2];
sx q[2];
rz(-0.85181063) q[2];
rz(-0.24438721) q[3];
sx q[3];
rz(-0.71072018) q[3];
sx q[3];
rz(-1.3062668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39563018) q[0];
sx q[0];
rz(-1.4816477) q[0];
sx q[0];
rz(-3.0356044) q[0];
rz(1.7164187) q[1];
sx q[1];
rz(-1.1499848) q[1];
sx q[1];
rz(-1.0894159) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46360923) q[0];
sx q[0];
rz(-2.8937723) q[0];
sx q[0];
rz(1.8576966) q[0];
rz(-pi) q[1];
rz(2.0136497) q[2];
sx q[2];
rz(-0.93292716) q[2];
sx q[2];
rz(1.4952212) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.7287561) q[1];
sx q[1];
rz(-2.7979971) q[1];
sx q[1];
rz(1.1854002) q[1];
rz(-1.7884729) q[3];
sx q[3];
rz(-2.5158415) q[3];
sx q[3];
rz(-1.1661287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0561698) q[2];
sx q[2];
rz(-1.9499754) q[2];
sx q[2];
rz(-2.5540111) q[2];
rz(-1.9116631) q[3];
sx q[3];
rz(-2.7432224) q[3];
sx q[3];
rz(0.57268322) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2040445) q[0];
sx q[0];
rz(-1.2296822) q[0];
sx q[0];
rz(-2.6958534) q[0];
rz(-0.84272376) q[1];
sx q[1];
rz(-2.3363967) q[1];
sx q[1];
rz(2.3420948) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3955376) q[0];
sx q[0];
rz(-2.5776401) q[0];
sx q[0];
rz(-0.82036316) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72321524) q[2];
sx q[2];
rz(-0.84645459) q[2];
sx q[2];
rz(-2.3220313) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.013019907) q[1];
sx q[1];
rz(-1.1718084) q[1];
sx q[1];
rz(0.47069957) q[1];
rz(-1.4700674) q[3];
sx q[3];
rz(-2.5593649) q[3];
sx q[3];
rz(1.738172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2004956) q[2];
sx q[2];
rz(-0.37797394) q[2];
sx q[2];
rz(-2.359158) q[2];
rz(-0.93067074) q[3];
sx q[3];
rz(-1.2717671) q[3];
sx q[3];
rz(-2.8455287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3928669) q[0];
sx q[0];
rz(-2.1935232) q[0];
sx q[0];
rz(-2.0178846) q[0];
rz(0.8404845) q[1];
sx q[1];
rz(-1.1437462) q[1];
sx q[1];
rz(-1.3927654) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3884149) q[0];
sx q[0];
rz(-1.6196189) q[0];
sx q[0];
rz(1.6724259) q[0];
rz(-pi) q[1];
rz(2.8181067) q[2];
sx q[2];
rz(-1.3703823) q[2];
sx q[2];
rz(-0.85497626) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1906793) q[1];
sx q[1];
rz(-1.9164521) q[1];
sx q[1];
rz(0.18195621) q[1];
x q[2];
rz(-1.6540307) q[3];
sx q[3];
rz(-2.1862767) q[3];
sx q[3];
rz(-2.7546101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.30354083) q[2];
sx q[2];
rz(-2.3736062) q[2];
sx q[2];
rz(1.6221907) q[2];
rz(1.0649902) q[3];
sx q[3];
rz(-2.0678554) q[3];
sx q[3];
rz(-0.85552335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72724718) q[0];
sx q[0];
rz(-2.1942744) q[0];
sx q[0];
rz(-1.974768) q[0];
rz(-2.836152) q[1];
sx q[1];
rz(-2.0716397) q[1];
sx q[1];
rz(0.50183141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0286642) q[0];
sx q[0];
rz(-1.477286) q[0];
sx q[0];
rz(-0.19569719) q[0];
x q[1];
rz(-0.86595777) q[2];
sx q[2];
rz(-2.188345) q[2];
sx q[2];
rz(2.7028529) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33602411) q[1];
sx q[1];
rz(-0.37266392) q[1];
sx q[1];
rz(-1.3678846) q[1];
x q[2];
rz(-2.3324729) q[3];
sx q[3];
rz(-0.71965766) q[3];
sx q[3];
rz(-1.9687195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66494232) q[2];
sx q[2];
rz(-2.3822337) q[2];
sx q[2];
rz(3.0090028) q[2];
rz(2.8442123) q[3];
sx q[3];
rz(-1.6011651) q[3];
sx q[3];
rz(2.4471349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7731758) q[0];
sx q[0];
rz(-1.6650124) q[0];
sx q[0];
rz(-1.6126527) q[0];
rz(-1.7932786) q[1];
sx q[1];
rz(-1.3765843) q[1];
sx q[1];
rz(-2.7682159) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7517029) q[0];
sx q[0];
rz(-1.9284964) q[0];
sx q[0];
rz(0.16660868) q[0];
x q[1];
rz(-0.53826314) q[2];
sx q[2];
rz(-0.54780947) q[2];
sx q[2];
rz(-2.4108568) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8638525) q[1];
sx q[1];
rz(-1.9406295) q[1];
sx q[1];
rz(-0.86994008) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67383234) q[3];
sx q[3];
rz(-0.34075173) q[3];
sx q[3];
rz(1.0125404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1030964) q[2];
sx q[2];
rz(-3.0615443) q[2];
sx q[2];
rz(2.6888964) q[2];
rz(-2.2321841) q[3];
sx q[3];
rz(-2.1287287) q[3];
sx q[3];
rz(0.67210853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4368923) q[0];
sx q[0];
rz(-1.1672651) q[0];
sx q[0];
rz(1.1080909) q[0];
rz(2.3469901) q[1];
sx q[1];
rz(-1.662685) q[1];
sx q[1];
rz(-0.86659238) q[1];
rz(1.5431326) q[2];
sx q[2];
rz(-2.4504708) q[2];
sx q[2];
rz(-2.8219555) q[2];
rz(-1.1232212) q[3];
sx q[3];
rz(-2.2991857) q[3];
sx q[3];
rz(-3.1274336) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
