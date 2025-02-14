OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95028967) q[0];
sx q[0];
rz(-0.27009717) q[0];
sx q[0];
rz(0.88859963) q[0];
rz(-1.3130045) q[1];
sx q[1];
rz(-1.5994025) q[1];
sx q[1];
rz(1.7607652) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.439405) q[0];
sx q[0];
rz(-2.8719465) q[0];
sx q[0];
rz(-2.3764971) q[0];
rz(-0.51989748) q[2];
sx q[2];
rz(-2.2714104) q[2];
sx q[2];
rz(-2.0127279) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1596229) q[1];
sx q[1];
rz(-1.2352984) q[1];
sx q[1];
rz(-1.580975) q[1];
x q[2];
rz(-2.9167261) q[3];
sx q[3];
rz(-1.2149842) q[3];
sx q[3];
rz(1.0575305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.31972739) q[2];
sx q[2];
rz(-1.3250019) q[2];
sx q[2];
rz(2.3925609) q[2];
rz(0.15394112) q[3];
sx q[3];
rz(-1.0356244) q[3];
sx q[3];
rz(3.0416987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6946436) q[0];
sx q[0];
rz(-1.6634989) q[0];
sx q[0];
rz(1.2167759) q[0];
rz(1.0617537) q[1];
sx q[1];
rz(-2.3371425) q[1];
sx q[1];
rz(0.42713508) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1841284) q[0];
sx q[0];
rz(-1.2250568) q[0];
sx q[0];
rz(2.5347802) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2160276) q[2];
sx q[2];
rz(-1.3750018) q[2];
sx q[2];
rz(-2.2044971) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4915575) q[1];
sx q[1];
rz(-1.2833529) q[1];
sx q[1];
rz(-0.14658714) q[1];
rz(-pi) q[2];
rz(0.052316523) q[3];
sx q[3];
rz(-1.3671095) q[3];
sx q[3];
rz(1.8546722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0624258) q[2];
sx q[2];
rz(-1.7288952) q[2];
sx q[2];
rz(-1.742935) q[2];
rz(2.9178197) q[3];
sx q[3];
rz(-2.2327435) q[3];
sx q[3];
rz(1.5435425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021521213) q[0];
sx q[0];
rz(-1.3503617) q[0];
sx q[0];
rz(0.045510005) q[0];
rz(1.1687763) q[1];
sx q[1];
rz(-1.903542) q[1];
sx q[1];
rz(-0.57560903) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5110014) q[0];
sx q[0];
rz(-1.5596034) q[0];
sx q[0];
rz(2.5024947) q[0];
x q[1];
rz(-0.91012886) q[2];
sx q[2];
rz(-1.5572539) q[2];
sx q[2];
rz(-2.8000268) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61316669) q[1];
sx q[1];
rz(-1.7878782) q[1];
sx q[1];
rz(-2.1291127) q[1];
rz(-pi) q[2];
rz(0.024954114) q[3];
sx q[3];
rz(-0.91812274) q[3];
sx q[3];
rz(-3.0252473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0098972926) q[2];
sx q[2];
rz(-0.74247777) q[2];
sx q[2];
rz(-2.1843145) q[2];
rz(-0.57473985) q[3];
sx q[3];
rz(-1.6114019) q[3];
sx q[3];
rz(1.8360229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9855373) q[0];
sx q[0];
rz(-2.1888581) q[0];
sx q[0];
rz(1.8804469) q[0];
rz(-0.082322923) q[1];
sx q[1];
rz(-2.0539093) q[1];
sx q[1];
rz(1.3538768) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.026431) q[0];
sx q[0];
rz(-1.8055834) q[0];
sx q[0];
rz(-2.874766) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33892314) q[2];
sx q[2];
rz(-0.51033516) q[2];
sx q[2];
rz(-1.5938544) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2980068) q[1];
sx q[1];
rz(-1.6248967) q[1];
sx q[1];
rz(-0.49523103) q[1];
rz(2.2500751) q[3];
sx q[3];
rz(-1.1783021) q[3];
sx q[3];
rz(-2.6160927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.42133078) q[2];
sx q[2];
rz(-0.91840863) q[2];
sx q[2];
rz(-1.8850231) q[2];
rz(0.71150696) q[3];
sx q[3];
rz(-1.810775) q[3];
sx q[3];
rz(2.8594657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8307777) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(-2.7984483) q[0];
rz(-0.061773069) q[1];
sx q[1];
rz(-0.97007483) q[1];
sx q[1];
rz(1.7247346) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8648711) q[0];
sx q[0];
rz(-1.7407932) q[0];
sx q[0];
rz(1.2890588) q[0];
x q[1];
rz(2.7602536) q[2];
sx q[2];
rz(-2.4040607) q[2];
sx q[2];
rz(-0.87070891) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7570863) q[1];
sx q[1];
rz(-2.050274) q[1];
sx q[1];
rz(0.15285413) q[1];
rz(-1.1457025) q[3];
sx q[3];
rz(-1.3561965) q[3];
sx q[3];
rz(-0.44073013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5938277) q[2];
sx q[2];
rz(-1.4676899) q[2];
sx q[2];
rz(2.055638) q[2];
rz(2.2222399) q[3];
sx q[3];
rz(-1.3708401) q[3];
sx q[3];
rz(-0.61387387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3980961) q[0];
sx q[0];
rz(-1.2092104) q[0];
sx q[0];
rz(0.79291517) q[0];
rz(-2.4010557) q[1];
sx q[1];
rz(-2.1426327) q[1];
sx q[1];
rz(-0.83121306) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8342469) q[0];
sx q[0];
rz(-1.9458658) q[0];
sx q[0];
rz(0.22386472) q[0];
rz(1.2937282) q[2];
sx q[2];
rz(-1.5779621) q[2];
sx q[2];
rz(-2.9286419) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.97396353) q[1];
sx q[1];
rz(-2.2669753) q[1];
sx q[1];
rz(0.40145282) q[1];
rz(1.5300203) q[3];
sx q[3];
rz(-0.47104657) q[3];
sx q[3];
rz(-0.99179964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.071216019) q[2];
sx q[2];
rz(-1.2717609) q[2];
sx q[2];
rz(2.0224723) q[2];
rz(-1.2707155) q[3];
sx q[3];
rz(-2.0344574) q[3];
sx q[3];
rz(1.3814111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1374461) q[0];
sx q[0];
rz(-2.9747712) q[0];
sx q[0];
rz(1.5420089) q[0];
rz(1.0150602) q[1];
sx q[1];
rz(-1.5559745) q[1];
sx q[1];
rz(2.9687845) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5046523) q[0];
sx q[0];
rz(-2.0604134) q[0];
sx q[0];
rz(2.2235653) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44266959) q[2];
sx q[2];
rz(-0.73390642) q[2];
sx q[2];
rz(-1.1709605) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1041278) q[1];
sx q[1];
rz(-1.8191254) q[1];
sx q[1];
rz(-2.009997) q[1];
x q[2];
rz(2.1234346) q[3];
sx q[3];
rz(-0.88200906) q[3];
sx q[3];
rz(1.3017201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66199866) q[2];
sx q[2];
rz(-2.2981503) q[2];
sx q[2];
rz(0.97770989) q[2];
rz(0.57502037) q[3];
sx q[3];
rz(-1.2049371) q[3];
sx q[3];
rz(1.9112816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8136895) q[0];
sx q[0];
rz(-2.8459025) q[0];
sx q[0];
rz(-3.1258702) q[0];
rz(-1.9533336) q[1];
sx q[1];
rz(-2.5091722) q[1];
sx q[1];
rz(1.1994919) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48391438) q[0];
sx q[0];
rz(-1.3356613) q[0];
sx q[0];
rz(-1.1170618) q[0];
rz(-pi) q[1];
rz(-2.0616777) q[2];
sx q[2];
rz(-1.617031) q[2];
sx q[2];
rz(-2.2904928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4839061) q[1];
sx q[1];
rz(-0.71345854) q[1];
sx q[1];
rz(-2.3785527) q[1];
rz(-2.0301129) q[3];
sx q[3];
rz(-2.5190341) q[3];
sx q[3];
rz(2.7186269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99379313) q[2];
sx q[2];
rz(-1.2387929) q[2];
sx q[2];
rz(-1.3851059) q[2];
rz(-0.6238474) q[3];
sx q[3];
rz(-1.2522937) q[3];
sx q[3];
rz(2.0417716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1856336) q[0];
sx q[0];
rz(-0.82056844) q[0];
sx q[0];
rz(-0.69325915) q[0];
rz(-0.26501003) q[1];
sx q[1];
rz(-2.3131504) q[1];
sx q[1];
rz(1.6185435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5361621) q[0];
sx q[0];
rz(-2.3443065) q[0];
sx q[0];
rz(1.5553586) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9829282) q[2];
sx q[2];
rz(-1.0318021) q[2];
sx q[2];
rz(-1.7807775) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0654782) q[1];
sx q[1];
rz(-1.5468742) q[1];
sx q[1];
rz(0.33965276) q[1];
rz(-3.070799) q[3];
sx q[3];
rz(-0.99081836) q[3];
sx q[3];
rz(-1.1759315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8388464) q[2];
sx q[2];
rz(-1.6281717) q[2];
sx q[2];
rz(1.6839074) q[2];
rz(0.95528209) q[3];
sx q[3];
rz(-0.49946076) q[3];
sx q[3];
rz(0.83520755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7569698) q[0];
sx q[0];
rz(-2.5769825) q[0];
sx q[0];
rz(-2.6589174) q[0];
rz(-0.92974281) q[1];
sx q[1];
rz(-1.0044121) q[1];
sx q[1];
rz(2.7453056) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1866465) q[0];
sx q[0];
rz(-0.56342331) q[0];
sx q[0];
rz(-2.4231829) q[0];
rz(-pi) q[1];
rz(1.4192788) q[2];
sx q[2];
rz(-1.3938402) q[2];
sx q[2];
rz(-2.3529476) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.63782802) q[1];
sx q[1];
rz(-2.2866268) q[1];
sx q[1];
rz(0.24055918) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9307053) q[3];
sx q[3];
rz(-1.1904089) q[3];
sx q[3];
rz(-1.4402657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79260176) q[2];
sx q[2];
rz(-1.7020117) q[2];
sx q[2];
rz(-0.3271884) q[2];
rz(0.78091019) q[3];
sx q[3];
rz(-1.1612929) q[3];
sx q[3];
rz(1.6646632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0384211) q[0];
sx q[0];
rz(-2.1506943) q[0];
sx q[0];
rz(-2.8839169) q[0];
rz(0.36956638) q[1];
sx q[1];
rz(-2.259544) q[1];
sx q[1];
rz(2.4824711) q[1];
rz(1.5255899) q[2];
sx q[2];
rz(-2.2635985) q[2];
sx q[2];
rz(0.23512693) q[2];
rz(1.0450324) q[3];
sx q[3];
rz(-0.67787328) q[3];
sx q[3];
rz(-1.9280435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
