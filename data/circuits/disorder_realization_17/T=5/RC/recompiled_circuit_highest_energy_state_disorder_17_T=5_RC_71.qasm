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
rz(2.0789335) q[0];
sx q[0];
rz(-2.3261429) q[0];
sx q[0];
rz(2.4315779) q[0];
rz(2.8275936) q[1];
sx q[1];
rz(-2.2065838) q[1];
sx q[1];
rz(-1.8097872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3478617) q[0];
sx q[0];
rz(-1.2958741) q[0];
sx q[0];
rz(1.8040276) q[0];
rz(-pi) q[1];
rz(2.7983886) q[2];
sx q[2];
rz(-1.3319974) q[2];
sx q[2];
rz(2.3607139) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1430681) q[1];
sx q[1];
rz(-1.3245163) q[1];
sx q[1];
rz(0.95264901) q[1];
x q[2];
rz(2.136257) q[3];
sx q[3];
rz(-2.4255803) q[3];
sx q[3];
rz(-0.93040066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.72267246) q[2];
sx q[2];
rz(-1.2166497) q[2];
sx q[2];
rz(-2.326272) q[2];
rz(-2.3319862) q[3];
sx q[3];
rz(-1.4987192) q[3];
sx q[3];
rz(2.0577551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0631436) q[0];
sx q[0];
rz(-2.0922631) q[0];
sx q[0];
rz(0.11058841) q[0];
rz(-0.94509205) q[1];
sx q[1];
rz(-2.7022305) q[1];
sx q[1];
rz(-1.5792712) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54540578) q[0];
sx q[0];
rz(-2.2431787) q[0];
sx q[0];
rz(-0.59242146) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1223824) q[2];
sx q[2];
rz(-1.9539991) q[2];
sx q[2];
rz(-2.2926083) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.64792127) q[1];
sx q[1];
rz(-0.4954547) q[1];
sx q[1];
rz(0.22101553) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0575297) q[3];
sx q[3];
rz(-0.81765122) q[3];
sx q[3];
rz(-0.14150894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9940146) q[2];
sx q[2];
rz(-0.14573228) q[2];
sx q[2];
rz(2.6914524) q[2];
rz(-2.0077997) q[3];
sx q[3];
rz(-1.6885875) q[3];
sx q[3];
rz(-2.1639737) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5728773) q[0];
sx q[0];
rz(-2.1935538) q[0];
sx q[0];
rz(0.28133389) q[0];
rz(1.7549134) q[1];
sx q[1];
rz(-1.4298871) q[1];
sx q[1];
rz(2.5414355) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48334405) q[0];
sx q[0];
rz(-1.6172101) q[0];
sx q[0];
rz(0.29057403) q[0];
x q[1];
rz(2.1776206) q[2];
sx q[2];
rz(-1.5247048) q[2];
sx q[2];
rz(2.5986236) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8474947) q[1];
sx q[1];
rz(-1.0007684) q[1];
sx q[1];
rz(-1.9600201) q[1];
x q[2];
rz(-0.5731606) q[3];
sx q[3];
rz(-1.7148682) q[3];
sx q[3];
rz(0.39502963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0448138) q[2];
sx q[2];
rz(-2.009095) q[2];
sx q[2];
rz(-2.0951994) q[2];
rz(0.3977631) q[3];
sx q[3];
rz(-2.4989276) q[3];
sx q[3];
rz(-3.0090295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1716877) q[0];
sx q[0];
rz(-3.1145018) q[0];
sx q[0];
rz(1.0525674) q[0];
rz(-1.196208) q[1];
sx q[1];
rz(-1.2095249) q[1];
sx q[1];
rz(2.722091) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9297816) q[0];
sx q[0];
rz(-0.63840862) q[0];
sx q[0];
rz(-2.6891461) q[0];
rz(-pi) q[1];
rz(-0.93650903) q[2];
sx q[2];
rz(-1.3872996) q[2];
sx q[2];
rz(-1.8734135) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4070523) q[1];
sx q[1];
rz(-0.95065763) q[1];
sx q[1];
rz(-2.6728515) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2798843) q[3];
sx q[3];
rz(-2.0613163) q[3];
sx q[3];
rz(1.377493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18998751) q[2];
sx q[2];
rz(-2.0230468) q[2];
sx q[2];
rz(-1.0901964) q[2];
rz(-0.53127855) q[3];
sx q[3];
rz(-0.71610206) q[3];
sx q[3];
rz(-1.6147015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4215609) q[0];
sx q[0];
rz(-1.5721385) q[0];
sx q[0];
rz(1.3979727) q[0];
rz(3.0086503) q[1];
sx q[1];
rz(-1.839919) q[1];
sx q[1];
rz(-3.1403819) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0069790445) q[0];
sx q[0];
rz(-1.190257) q[0];
sx q[0];
rz(-1.5534205) q[0];
x q[1];
rz(-2.9820061) q[2];
sx q[2];
rz(-1.8529466) q[2];
sx q[2];
rz(0.42195502) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7442786) q[1];
sx q[1];
rz(-0.71956735) q[1];
sx q[1];
rz(2.6707763) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7453543) q[3];
sx q[3];
rz(-2.9896185) q[3];
sx q[3];
rz(2.0106955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.91263897) q[2];
sx q[2];
rz(-1.3372083) q[2];
sx q[2];
rz(-0.21052989) q[2];
rz(1.0362961) q[3];
sx q[3];
rz(-0.784289) q[3];
sx q[3];
rz(-1.2150631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1933111) q[0];
sx q[0];
rz(-1.9183777) q[0];
sx q[0];
rz(2.7358828) q[0];
rz(-1.3544719) q[1];
sx q[1];
rz(-0.82258075) q[1];
sx q[1];
rz(-2.2192661) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2012537) q[0];
sx q[0];
rz(-2.8329599) q[0];
sx q[0];
rz(-1.2135394) q[0];
rz(-0.90409235) q[2];
sx q[2];
rz(-1.4897926) q[2];
sx q[2];
rz(1.9503649) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4757918) q[1];
sx q[1];
rz(-2.5195055) q[1];
sx q[1];
rz(1.3740963) q[1];
rz(-pi) q[2];
rz(-2.5821463) q[3];
sx q[3];
rz(-1.2370951) q[3];
sx q[3];
rz(-2.4041686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4313844) q[2];
sx q[2];
rz(-1.9269383) q[2];
sx q[2];
rz(0.40194884) q[2];
rz(-1.8534144) q[3];
sx q[3];
rz(-0.86789075) q[3];
sx q[3];
rz(-1.8624381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61889082) q[0];
sx q[0];
rz(-2.0382477) q[0];
sx q[0];
rz(3.0772305) q[0];
rz(-2.3035658) q[1];
sx q[1];
rz(-2.3152654) q[1];
sx q[1];
rz(-1.3305957) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15312402) q[0];
sx q[0];
rz(-2.5543384) q[0];
sx q[0];
rz(-2.9140317) q[0];
rz(-1.072791) q[2];
sx q[2];
rz(-2.118022) q[2];
sx q[2];
rz(-2.4233873) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8992675) q[1];
sx q[1];
rz(-1.0501672) q[1];
sx q[1];
rz(-1.316687) q[1];
rz(0.17223151) q[3];
sx q[3];
rz(-2.1759998) q[3];
sx q[3];
rz(1.8021405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4055206) q[2];
sx q[2];
rz(-1.5084927) q[2];
sx q[2];
rz(2.4580477) q[2];
rz(-2.4077967) q[3];
sx q[3];
rz(-1.7283231) q[3];
sx q[3];
rz(2.2939513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7965294) q[0];
sx q[0];
rz(-2.0579484) q[0];
sx q[0];
rz(-1.1731359) q[0];
rz(0.37551156) q[1];
sx q[1];
rz(-0.48986062) q[1];
sx q[1];
rz(-1.0367905) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4597367) q[0];
sx q[0];
rz(-1.5516743) q[0];
sx q[0];
rz(-0.023650344) q[0];
rz(-pi) q[1];
rz(-2.9762245) q[2];
sx q[2];
rz(-1.2980808) q[2];
sx q[2];
rz(-0.80406666) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.45932508) q[1];
sx q[1];
rz(-0.81574355) q[1];
sx q[1];
rz(1.3774648) q[1];
x q[2];
rz(-2.6508207) q[3];
sx q[3];
rz(-2.0476855) q[3];
sx q[3];
rz(-2.8618438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.57572395) q[2];
sx q[2];
rz(-0.91369358) q[2];
sx q[2];
rz(0.73053989) q[2];
rz(0.20914397) q[3];
sx q[3];
rz(-2.6181965) q[3];
sx q[3];
rz(1.8714347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4850979) q[0];
sx q[0];
rz(-0.42335835) q[0];
sx q[0];
rz(3.094161) q[0];
rz(0.221953) q[1];
sx q[1];
rz(-1.0739948) q[1];
sx q[1];
rz(-0.91112959) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.52235) q[0];
sx q[0];
rz(-1.5372835) q[0];
sx q[0];
rz(-1.5438118) q[0];
x q[1];
rz(-2.2124771) q[2];
sx q[2];
rz(-1.5589899) q[2];
sx q[2];
rz(-0.75071834) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4998656) q[1];
sx q[1];
rz(-1.9382825) q[1];
sx q[1];
rz(-1.5572474) q[1];
rz(-pi) q[2];
rz(1.2558612) q[3];
sx q[3];
rz(-0.58940998) q[3];
sx q[3];
rz(-2.3737597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.86965108) q[2];
sx q[2];
rz(-0.38001529) q[2];
sx q[2];
rz(0.16858777) q[2];
rz(2.9224959) q[3];
sx q[3];
rz(-0.85182652) q[3];
sx q[3];
rz(1.3271837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.99437) q[0];
sx q[0];
rz(-2.8012025) q[0];
sx q[0];
rz(0.45743531) q[0];
rz(-1.9153197) q[1];
sx q[1];
rz(-0.33718449) q[1];
sx q[1];
rz(-1.7785243) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5209893) q[0];
sx q[0];
rz(-2.531157) q[0];
sx q[0];
rz(0.10535289) q[0];
x q[1];
rz(0.3293475) q[2];
sx q[2];
rz(-2.4852356) q[2];
sx q[2];
rz(-2.658297) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.05039617) q[1];
sx q[1];
rz(-1.947787) q[1];
sx q[1];
rz(-1.7524377) q[1];
x q[2];
rz(-0.8882723) q[3];
sx q[3];
rz(-1.5049045) q[3];
sx q[3];
rz(-0.94277387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7676131) q[2];
sx q[2];
rz(-2.5869936) q[2];
sx q[2];
rz(0.45200959) q[2];
rz(-0.40397817) q[3];
sx q[3];
rz(-1.890506) q[3];
sx q[3];
rz(2.0154791) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44611888) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(-0.52234621) q[1];
sx q[1];
rz(-0.53032395) q[1];
sx q[1];
rz(1.2293336) q[1];
rz(1.8277373) q[2];
sx q[2];
rz(-1.8443454) q[2];
sx q[2];
rz(0.92583427) q[2];
rz(1.4215076) q[3];
sx q[3];
rz(-1.310077) q[3];
sx q[3];
rz(2.6826912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
