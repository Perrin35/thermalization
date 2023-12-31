OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.26602715) q[0];
sx q[0];
rz(-0.53524435) q[0];
sx q[0];
rz(0.75403655) q[0];
rz(0.79020483) q[1];
sx q[1];
rz(-1.2269998) q[1];
sx q[1];
rz(-1.1608646) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0275118) q[0];
sx q[0];
rz(-2.049276) q[0];
sx q[0];
rz(1.6041243) q[0];
rz(-pi) q[1];
rz(1.8765175) q[2];
sx q[2];
rz(-2.6439878) q[2];
sx q[2];
rz(1.9203609) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.3996934) q[1];
sx q[1];
rz(-1.0642991) q[1];
sx q[1];
rz(-2.6687117) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74673064) q[3];
sx q[3];
rz(-0.24198469) q[3];
sx q[3];
rz(-0.63499588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6291818) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(-0.66649246) q[2];
rz(-0.50981057) q[3];
sx q[3];
rz(-1.1504268) q[3];
sx q[3];
rz(1.2734909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3614685) q[0];
sx q[0];
rz(-1.7391917) q[0];
sx q[0];
rz(2.0289039) q[0];
rz(-0.15377046) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(1.8033093) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7897658) q[0];
sx q[0];
rz(-1.7625945) q[0];
sx q[0];
rz(0.70621323) q[0];
x q[1];
rz(2.0427225) q[2];
sx q[2];
rz(-2.3510691) q[2];
sx q[2];
rz(-2.3015442) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2087304) q[1];
sx q[1];
rz(-1.6881697) q[1];
sx q[1];
rz(-2.7670303) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5541833) q[3];
sx q[3];
rz(-1.4360768) q[3];
sx q[3];
rz(-0.87232529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0587557) q[2];
sx q[2];
rz(-2.357491) q[2];
sx q[2];
rz(1.9632957) q[2];
rz(0.96238771) q[3];
sx q[3];
rz(-1.0932086) q[3];
sx q[3];
rz(0.66550955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798379) q[0];
sx q[0];
rz(-0.50757718) q[0];
sx q[0];
rz(1.6145561) q[0];
rz(-0.64287341) q[1];
sx q[1];
rz(-2.0650654) q[1];
sx q[1];
rz(0.33338526) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9133559) q[0];
sx q[0];
rz(-2.2524999) q[0];
sx q[0];
rz(2.8744254) q[0];
rz(1.2522167) q[2];
sx q[2];
rz(-0.40553667) q[2];
sx q[2];
rz(-0.47682724) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4585321) q[1];
sx q[1];
rz(-1.6525869) q[1];
sx q[1];
rz(-1.0973147) q[1];
rz(-pi) q[2];
rz(-0.63342996) q[3];
sx q[3];
rz(-1.389635) q[3];
sx q[3];
rz(1.0599979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12144111) q[2];
sx q[2];
rz(-1.3139498) q[2];
sx q[2];
rz(-0.80231673) q[2];
rz(2.9004167) q[3];
sx q[3];
rz(-2.4427588) q[3];
sx q[3];
rz(-3.0407217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0543095) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(2.5774082) q[0];
rz(0.57812771) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(0.50813466) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0429338) q[0];
sx q[0];
rz(-1.9693976) q[0];
sx q[0];
rz(-1.1111141) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5925243) q[2];
sx q[2];
rz(-0.93291514) q[2];
sx q[2];
rz(-2.4456172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9307738) q[1];
sx q[1];
rz(-0.79508077) q[1];
sx q[1];
rz(0.37097431) q[1];
x q[2];
rz(-1.2933613) q[3];
sx q[3];
rz(-1.7460896) q[3];
sx q[3];
rz(-0.41269916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6886787) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-2.508146) q[2];
rz(2.5417035) q[3];
sx q[3];
rz(-1.9918631) q[3];
sx q[3];
rz(1.6413123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7552345) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(0.19609837) q[0];
rz(1.261699) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(2.0702147) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8349985) q[0];
sx q[0];
rz(-2.1822565) q[0];
sx q[0];
rz(0.5284662) q[0];
x q[1];
rz(-1.1159665) q[2];
sx q[2];
rz(-2.358846) q[2];
sx q[2];
rz(0.46403971) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2744659) q[1];
sx q[1];
rz(-1.9879838) q[1];
sx q[1];
rz(-0.31068128) q[1];
rz(2.7171633) q[3];
sx q[3];
rz(-1.9924581) q[3];
sx q[3];
rz(-1.9205586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1065958) q[2];
sx q[2];
rz(-1.3369766) q[2];
sx q[2];
rz(2.1595947) q[2];
rz(-0.18520959) q[3];
sx q[3];
rz(-2.2976112) q[3];
sx q[3];
rz(1.265032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.5380602) q[0];
sx q[0];
rz(-2.2221727) q[0];
sx q[0];
rz(0.53034267) q[0];
rz(1.2999339) q[1];
sx q[1];
rz(-1.8118186) q[1];
sx q[1];
rz(-2.9249654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075561698) q[0];
sx q[0];
rz(-0.78853411) q[0];
sx q[0];
rz(1.1939474) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5112721) q[2];
sx q[2];
rz(-1.8152871) q[2];
sx q[2];
rz(-0.47843168) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4649012) q[1];
sx q[1];
rz(-1.5790107) q[1];
sx q[1];
rz(-1.2354922) q[1];
rz(-0.90029193) q[3];
sx q[3];
rz(-1.4196463) q[3];
sx q[3];
rz(-0.33526648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.65016046) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6525456) q[0];
sx q[0];
rz(-1.9788195) q[0];
sx q[0];
rz(0.52893692) q[0];
rz(-1.6128929) q[1];
sx q[1];
rz(-1.1922319) q[1];
sx q[1];
rz(1.0891917) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3952336) q[0];
sx q[0];
rz(-2.8440209) q[0];
sx q[0];
rz(1.5514821) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5838926) q[2];
sx q[2];
rz(-0.93765646) q[2];
sx q[2];
rz(1.5483088) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4186801) q[1];
sx q[1];
rz(-1.5376248) q[1];
sx q[1];
rz(1.5772217) q[1];
rz(-3.0508556) q[3];
sx q[3];
rz(-1.4416579) q[3];
sx q[3];
rz(-1.7307626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0163394) q[2];
sx q[2];
rz(-2.3374127) q[2];
sx q[2];
rz(-1.2255229) q[2];
rz(-1.6493753) q[3];
sx q[3];
rz(-1.4368613) q[3];
sx q[3];
rz(2.9857181) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65687031) q[0];
sx q[0];
rz(-3.0796034) q[0];
sx q[0];
rz(2.2739676) q[0];
rz(3.0743657) q[1];
sx q[1];
rz(-2.1104689) q[1];
sx q[1];
rz(0.19518383) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.09571) q[0];
sx q[0];
rz(-1.7716496) q[0];
sx q[0];
rz(0.30028371) q[0];
rz(2.7404501) q[2];
sx q[2];
rz(-1.5578798) q[2];
sx q[2];
rz(2.0502979) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.47123779) q[1];
sx q[1];
rz(-2.0002623) q[1];
sx q[1];
rz(-3.0958789) q[1];
rz(-pi) q[2];
rz(-0.35555367) q[3];
sx q[3];
rz(-2.0204633) q[3];
sx q[3];
rz(-0.21330968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8937257) q[2];
sx q[2];
rz(-2.9276431) q[2];
sx q[2];
rz(0.90551886) q[2];
rz(1.1577822) q[3];
sx q[3];
rz(-1.1685305) q[3];
sx q[3];
rz(2.982443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72682196) q[0];
sx q[0];
rz(-1.0960217) q[0];
sx q[0];
rz(-2.1006405) q[0];
rz(-3.0629311) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(0.35531607) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2040877) q[0];
sx q[0];
rz(-1.4728439) q[0];
sx q[0];
rz(1.9393001) q[0];
rz(-pi) q[1];
rz(1.0742498) q[2];
sx q[2];
rz(-1.7830007) q[2];
sx q[2];
rz(-0.70763904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73341093) q[1];
sx q[1];
rz(-1.8942041) q[1];
sx q[1];
rz(3.0820993) q[1];
rz(-1.7570417) q[3];
sx q[3];
rz(-2.1631141) q[3];
sx q[3];
rz(0.67519644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67655247) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(-1.4036277) q[2];
rz(-0.0062395652) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(1.5374373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86826098) q[0];
sx q[0];
rz(-1.9737759) q[0];
sx q[0];
rz(-1.8433174) q[0];
rz(0.4459933) q[1];
sx q[1];
rz(-2.0263717) q[1];
sx q[1];
rz(0.30074063) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8200127) q[0];
sx q[0];
rz(-2.1029148) q[0];
sx q[0];
rz(2.8768455) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54603521) q[2];
sx q[2];
rz(-1.2166426) q[2];
sx q[2];
rz(1.491577) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.7584596) q[1];
sx q[1];
rz(-2.8499095) q[1];
sx q[1];
rz(-0.056706927) q[1];
x q[2];
rz(-2.2311287) q[3];
sx q[3];
rz(-1.904389) q[3];
sx q[3];
rz(3.0433082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62844244) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(-2.002031) q[2];
rz(-1.7906174) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(-1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.448485) q[0];
sx q[0];
rz(-1.9036475) q[0];
sx q[0];
rz(0.87686476) q[0];
rz(1.4032455) q[1];
sx q[1];
rz(-1.9156024) q[1];
sx q[1];
rz(1.8617873) q[1];
rz(0.055988978) q[2];
sx q[2];
rz(-1.5691139) q[2];
sx q[2];
rz(-0.72473095) q[2];
rz(2.8786447) q[3];
sx q[3];
rz(-2.1750952) q[3];
sx q[3];
rz(1.7261214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
