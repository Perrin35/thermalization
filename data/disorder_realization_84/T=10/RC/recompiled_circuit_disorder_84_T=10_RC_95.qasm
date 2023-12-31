OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8755655) q[0];
sx q[0];
rz(-2.6063483) q[0];
sx q[0];
rz(2.3875561) q[0];
rz(-5.4929805) q[1];
sx q[1];
rz(5.0561855) q[1];
sx q[1];
rz(8.2639134) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1140808) q[0];
sx q[0];
rz(-1.0923166) q[0];
sx q[0];
rz(1.6041243) q[0];
rz(-pi) q[1];
rz(-1.0928458) q[2];
sx q[2];
rz(-1.7149601) q[2];
sx q[2];
rz(0.62010566) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7293284) q[1];
sx q[1];
rz(-0.67854478) q[1];
sx q[1];
rz(-2.2581984) q[1];
x q[2];
rz(1.7369032) q[3];
sx q[3];
rz(-1.7475834) q[3];
sx q[3];
rz(-0.12648957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51241088) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(0.66649246) q[2];
rz(-2.6317821) q[3];
sx q[3];
rz(-1.9911659) q[3];
sx q[3];
rz(-1.8681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3614685) q[0];
sx q[0];
rz(-1.402401) q[0];
sx q[0];
rz(1.1126888) q[0];
rz(0.15377046) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(1.3382834) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99909821) q[0];
sx q[0];
rz(-0.72744838) q[0];
sx q[0];
rz(2.8508458) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0427225) q[2];
sx q[2];
rz(-0.79052351) q[2];
sx q[2];
rz(0.84004842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8255303) q[1];
sx q[1];
rz(-1.1989374) q[1];
sx q[1];
rz(-1.6968246) q[1];
x q[2];
rz(0.23985858) q[3];
sx q[3];
rz(-0.60088241) q[3];
sx q[3];
rz(0.89752737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0587557) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(-1.9632957) q[2];
rz(0.96238771) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(-0.66550955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798379) q[0];
sx q[0];
rz(-0.50757718) q[0];
sx q[0];
rz(-1.6145561) q[0];
rz(-0.64287341) q[1];
sx q[1];
rz(-2.0650654) q[1];
sx q[1];
rz(0.33338526) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3231455) q[0];
sx q[0];
rz(-0.72431699) q[0];
sx q[0];
rz(-1.2562654) q[0];
rz(1.2522167) q[2];
sx q[2];
rz(-0.40553667) q[2];
sx q[2];
rz(2.6647654) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4585321) q[1];
sx q[1];
rz(-1.4890057) q[1];
sx q[1];
rz(1.0973147) q[1];
rz(-0.63342996) q[3];
sx q[3];
rz(-1.7519577) q[3];
sx q[3];
rz(2.0815947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12144111) q[2];
sx q[2];
rz(-1.3139498) q[2];
sx q[2];
rz(0.80231673) q[2];
rz(-2.9004167) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(-3.0407217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0543095) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(-0.56418443) q[0];
rz(0.57812771) q[1];
sx q[1];
rz(-1.6491978) q[1];
sx q[1];
rz(-0.50813466) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098658853) q[0];
sx q[0];
rz(-1.172195) q[0];
sx q[0];
rz(2.0304785) q[0];
rz(-0.029301734) q[2];
sx q[2];
rz(-0.63819956) q[2];
sx q[2];
rz(-2.4820941) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9307738) q[1];
sx q[1];
rz(-0.79508077) q[1];
sx q[1];
rz(0.37097431) q[1];
rz(-pi) q[2];
rz(0.99677892) q[3];
sx q[3];
rz(-2.814631) q[3];
sx q[3];
rz(2.5329778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6886787) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(0.63344669) q[2];
rz(-2.5417035) q[3];
sx q[3];
rz(-1.9918631) q[3];
sx q[3];
rz(1.5002804) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7552345) q[0];
sx q[0];
rz(-2.8383377) q[0];
sx q[0];
rz(2.9454943) q[0];
rz(1.261699) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(2.0702147) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6288554) q[0];
sx q[0];
rz(-0.78539408) q[0];
sx q[0];
rz(0.94731713) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1159665) q[2];
sx q[2];
rz(-2.358846) q[2];
sx q[2];
rz(-0.46403971) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5672946) q[1];
sx q[1];
rz(-1.8540566) q[1];
sx q[1];
rz(1.1351372) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42442935) q[3];
sx q[3];
rz(-1.1491346) q[3];
sx q[3];
rz(-1.9205586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.034996899) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(-2.1595947) q[2];
rz(2.9563831) q[3];
sx q[3];
rz(-2.2976112) q[3];
sx q[3];
rz(1.265032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6035325) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-2.61125) q[0];
rz(-1.8416587) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(-0.21662724) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075561698) q[0];
sx q[0];
rz(-2.3530585) q[0];
sx q[0];
rz(-1.9476452) q[0];
rz(-pi) q[1];
rz(1.8703307) q[2];
sx q[2];
rz(-2.1795863) q[2];
sx q[2];
rz(-1.2671721) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12946373) q[1];
sx q[1];
rz(-2.8061917) q[1];
sx q[1];
rz(1.5458376) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8111749) q[3];
sx q[3];
rz(-0.68475311) q[3];
sx q[3];
rz(-1.7184337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4914322) q[2];
sx q[2];
rz(-1.5101134) q[2];
sx q[2];
rz(-2.1177297) q[2];
rz(-2.2653545) q[3];
sx q[3];
rz(-0.70228464) q[3];
sx q[3];
rz(1.2715626) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4890471) q[0];
sx q[0];
rz(-1.9788195) q[0];
sx q[0];
rz(-2.6126557) q[0];
rz(1.6128929) q[1];
sx q[1];
rz(-1.1922319) q[1];
sx q[1];
rz(-1.0891917) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.335621) q[0];
sx q[0];
rz(-1.5651337) q[0];
sx q[0];
rz(-1.8683158) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5577001) q[2];
sx q[2];
rz(-0.93765646) q[2];
sx q[2];
rz(1.5483088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4186801) q[1];
sx q[1];
rz(-1.6039679) q[1];
sx q[1];
rz(1.5772217) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7004622) q[3];
sx q[3];
rz(-1.6607758) q[3];
sx q[3];
rz(-2.9699096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0163394) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(1.2255229) q[2];
rz(1.6493753) q[3];
sx q[3];
rz(-1.4368613) q[3];
sx q[3];
rz(0.15587458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65687031) q[0];
sx q[0];
rz(-3.0796034) q[0];
sx q[0];
rz(0.86762506) q[0];
rz(3.0743657) q[1];
sx q[1];
rz(-2.1104689) q[1];
sx q[1];
rz(-2.9464088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0474437) q[0];
sx q[0];
rz(-0.35956811) q[0];
sx q[0];
rz(-2.5387499) q[0];
x q[1];
rz(-0.40114258) q[2];
sx q[2];
rz(-1.5578798) q[2];
sx q[2];
rz(-1.0912947) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47123779) q[1];
sx q[1];
rz(-2.0002623) q[1];
sx q[1];
rz(3.0958789) q[1];
x q[2];
rz(-2.1956452) q[3];
sx q[3];
rz(-2.5759856) q[3];
sx q[3];
rz(-2.2212976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.247867) q[2];
sx q[2];
rz(-2.9276431) q[2];
sx q[2];
rz(0.90551886) q[2];
rz(-1.1577822) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(-0.15914966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.72682196) q[0];
sx q[0];
rz(-1.0960217) q[0];
sx q[0];
rz(-2.1006405) q[0];
rz(3.0629311) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(2.7862766) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61475635) q[0];
sx q[0];
rz(-2.7608681) q[0];
sx q[0];
rz(1.8371131) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0742498) q[2];
sx q[2];
rz(-1.7830007) q[2];
sx q[2];
rz(0.70763904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3231343) q[1];
sx q[1];
rz(-1.5143906) q[1];
sx q[1];
rz(1.8947381) q[1];
rz(-pi) q[2];
rz(-2.5411685) q[3];
sx q[3];
rz(-1.416559) q[3];
sx q[3];
rz(-2.141181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.67655247) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(-1.4036277) q[2];
rz(0.0062395652) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(-1.5374373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733317) q[0];
sx q[0];
rz(-1.1678168) q[0];
sx q[0];
rz(-1.8433174) q[0];
rz(0.4459933) q[1];
sx q[1];
rz(-1.1152209) q[1];
sx q[1];
rz(-0.30074063) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3215799) q[0];
sx q[0];
rz(-1.0386779) q[0];
sx q[0];
rz(-0.26474712) q[0];
rz(-pi) q[1];
rz(-1.1624596) q[2];
sx q[2];
rz(-1.0620585) q[2];
sx q[2];
rz(-3.0131154) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69925752) q[1];
sx q[1];
rz(-1.8619969) q[1];
sx q[1];
rz(1.5878116) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91046393) q[3];
sx q[3];
rz(-1.2372036) q[3];
sx q[3];
rz(-3.0433082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5131502) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(-1.1395617) q[2];
rz(-1.3509753) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(-1.9679507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.448485) q[0];
sx q[0];
rz(-1.9036475) q[0];
sx q[0];
rz(0.87686476) q[0];
rz(1.7383472) q[1];
sx q[1];
rz(-1.2259903) q[1];
sx q[1];
rz(-1.2798053) q[1];
rz(1.5724814) q[2];
sx q[2];
rz(-1.6267852) q[2];
sx q[2];
rz(0.84597107) q[2];
rz(0.26294796) q[3];
sx q[3];
rz(-0.96649747) q[3];
sx q[3];
rz(-1.4154712) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
