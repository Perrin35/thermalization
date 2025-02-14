OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0176528) q[0];
sx q[0];
rz(4.058429) q[0];
sx q[0];
rz(9.8596758) q[0];
rz(1.5974367) q[1];
sx q[1];
rz(-2.6225852) q[1];
sx q[1];
rz(2.5595698) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30729025) q[0];
sx q[0];
rz(-0.30481642) q[0];
sx q[0];
rz(1.0340263) q[0];
x q[1];
rz(0.97442128) q[2];
sx q[2];
rz(-2.4108464) q[2];
sx q[2];
rz(1.7898909) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.21596244) q[1];
sx q[1];
rz(-2.0019889) q[1];
sx q[1];
rz(-0.27121065) q[1];
rz(-2.1277027) q[3];
sx q[3];
rz(-0.32809533) q[3];
sx q[3];
rz(-1.8522287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.55328289) q[2];
sx q[2];
rz(-1.8841691) q[2];
sx q[2];
rz(-0.29169875) q[2];
rz(2.6907673) q[3];
sx q[3];
rz(-2.8901633) q[3];
sx q[3];
rz(-2.5341471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70158231) q[0];
sx q[0];
rz(-2.1511183) q[0];
sx q[0];
rz(-1.095358) q[0];
rz(1.1913242) q[1];
sx q[1];
rz(-2.2215863) q[1];
sx q[1];
rz(-2.2390168) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1987225) q[0];
sx q[0];
rz(-2.1784221) q[0];
sx q[0];
rz(-0.81815079) q[0];
x q[1];
rz(1.9368725) q[2];
sx q[2];
rz(-0.70209018) q[2];
sx q[2];
rz(0.13167189) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.892977) q[1];
sx q[1];
rz(-1.6910403) q[1];
sx q[1];
rz(-2.8477232) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91530494) q[3];
sx q[3];
rz(-2.4065131) q[3];
sx q[3];
rz(2.8171768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9537182) q[2];
sx q[2];
rz(-0.97430054) q[2];
sx q[2];
rz(3.0226959) q[2];
rz(0.64905727) q[3];
sx q[3];
rz(-0.18638149) q[3];
sx q[3];
rz(-1.6868748) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4897937) q[0];
sx q[0];
rz(-1.2540023) q[0];
sx q[0];
rz(-2.6258262) q[0];
rz(2.0491397) q[1];
sx q[1];
rz(-2.7383883) q[1];
sx q[1];
rz(1.2752424) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3608801) q[0];
sx q[0];
rz(-3.096252) q[0];
sx q[0];
rz(-0.41098292) q[0];
rz(0.62531535) q[2];
sx q[2];
rz(-0.38039243) q[2];
sx q[2];
rz(1.8521295) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66652775) q[1];
sx q[1];
rz(-1.766886) q[1];
sx q[1];
rz(2.2580819) q[1];
x q[2];
rz(-1.6728064) q[3];
sx q[3];
rz(-2.3622741) q[3];
sx q[3];
rz(1.23097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0755997) q[2];
sx q[2];
rz(-1.3174572) q[2];
sx q[2];
rz(-1.9817748) q[2];
rz(2.6521111) q[3];
sx q[3];
rz(-1.5533841) q[3];
sx q[3];
rz(-0.71973962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8664261) q[0];
sx q[0];
rz(-0.32298276) q[0];
sx q[0];
rz(-0.48165709) q[0];
rz(0.9306759) q[1];
sx q[1];
rz(-1.8447256) q[1];
sx q[1];
rz(-3.1226588) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1279432) q[0];
sx q[0];
rz(-0.054220323) q[0];
sx q[0];
rz(-1.9677866) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2332814) q[2];
sx q[2];
rz(-1.5358972) q[2];
sx q[2];
rz(0.95373431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.035012951) q[1];
sx q[1];
rz(-1.0570608) q[1];
sx q[1];
rz(2.4519743) q[1];
rz(2.3778524) q[3];
sx q[3];
rz(-0.9943878) q[3];
sx q[3];
rz(1.1820205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21861741) q[2];
sx q[2];
rz(-0.3564035) q[2];
sx q[2];
rz(0.74529988) q[2];
rz(-0.37799147) q[3];
sx q[3];
rz(-1.9972921) q[3];
sx q[3];
rz(-0.88855827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9208263) q[0];
sx q[0];
rz(-2.8296318) q[0];
sx q[0];
rz(1.1908603) q[0];
rz(2.6530755) q[1];
sx q[1];
rz(-0.91594511) q[1];
sx q[1];
rz(0.8078422) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90276001) q[0];
sx q[0];
rz(-1.2441946) q[0];
sx q[0];
rz(0.095681698) q[0];
rz(-pi) q[1];
rz(2.9132782) q[2];
sx q[2];
rz(-2.0234152) q[2];
sx q[2];
rz(1.7646947) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3563429) q[1];
sx q[1];
rz(-1.6587573) q[1];
sx q[1];
rz(2.6846519) q[1];
rz(-2.1937815) q[3];
sx q[3];
rz(-1.4075507) q[3];
sx q[3];
rz(0.80006525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0025582) q[2];
sx q[2];
rz(-2.148874) q[2];
sx q[2];
rz(-0.16415088) q[2];
rz(0.74603355) q[3];
sx q[3];
rz(-1.371871) q[3];
sx q[3];
rz(0.073808864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1604851) q[0];
sx q[0];
rz(-1.8758513) q[0];
sx q[0];
rz(-0.72738457) q[0];
rz(-2.6630867) q[1];
sx q[1];
rz(-1.6886657) q[1];
sx q[1];
rz(-2.4047638) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7751303) q[0];
sx q[0];
rz(-1.7084165) q[0];
sx q[0];
rz(0.055995106) q[0];
x q[1];
rz(-0.13104266) q[2];
sx q[2];
rz(-1.2918345) q[2];
sx q[2];
rz(-2.758858) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2347654) q[1];
sx q[1];
rz(-1.2052457) q[1];
sx q[1];
rz(0.91616456) q[1];
rz(-pi) q[2];
rz(0.013929587) q[3];
sx q[3];
rz(-1.1831814) q[3];
sx q[3];
rz(2.882021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9937146) q[2];
sx q[2];
rz(-2.1174049) q[2];
sx q[2];
rz(2.4564339) q[2];
rz(-2.1327175) q[3];
sx q[3];
rz(-0.69005552) q[3];
sx q[3];
rz(-0.25079295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20550263) q[0];
sx q[0];
rz(-0.093955366) q[0];
sx q[0];
rz(2.0482546) q[0];
rz(-0.023748485) q[1];
sx q[1];
rz(-0.7370342) q[1];
sx q[1];
rz(-0.49560961) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8596397) q[0];
sx q[0];
rz(-3.0744327) q[0];
sx q[0];
rz(-0.18224506) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7772113) q[2];
sx q[2];
rz(-1.141618) q[2];
sx q[2];
rz(-0.51233722) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88265282) q[1];
sx q[1];
rz(-1.3892838) q[1];
sx q[1];
rz(-2.7869999) q[1];
x q[2];
rz(2.7109259) q[3];
sx q[3];
rz(-1.1209295) q[3];
sx q[3];
rz(-1.7783742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2451943) q[2];
sx q[2];
rz(-2.3485025) q[2];
sx q[2];
rz(-0.29655656) q[2];
rz(-0.5101997) q[3];
sx q[3];
rz(-1.5254131) q[3];
sx q[3];
rz(2.8701674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9940014) q[0];
sx q[0];
rz(-2.1864102) q[0];
sx q[0];
rz(-1.9872794) q[0];
rz(-1.1555903) q[1];
sx q[1];
rz(-2.4842333) q[1];
sx q[1];
rz(1.6824228) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4302368) q[0];
sx q[0];
rz(-0.85604307) q[0];
sx q[0];
rz(1.6011333) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1404678) q[2];
sx q[2];
rz(-0.84648593) q[2];
sx q[2];
rz(-1.520919) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1016804) q[1];
sx q[1];
rz(-1.3915359) q[1];
sx q[1];
rz(-2.0800566) q[1];
rz(-2.4637163) q[3];
sx q[3];
rz(-0.75784412) q[3];
sx q[3];
rz(2.0539396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3079754) q[2];
sx q[2];
rz(-2.5994382) q[2];
sx q[2];
rz(-1.3758434) q[2];
rz(3.0552676) q[3];
sx q[3];
rz(-0.83266801) q[3];
sx q[3];
rz(-1.8858006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1983222) q[0];
sx q[0];
rz(-0.69381303) q[0];
sx q[0];
rz(-0.86945239) q[0];
rz(2.4122639) q[1];
sx q[1];
rz(-0.84356934) q[1];
sx q[1];
rz(-2.9076911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6399891) q[0];
sx q[0];
rz(-2.7219048) q[0];
sx q[0];
rz(0.11735015) q[0];
rz(0.21017615) q[2];
sx q[2];
rz(-1.6646204) q[2];
sx q[2];
rz(-1.891234) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6944551) q[1];
sx q[1];
rz(-1.7181953) q[1];
sx q[1];
rz(-0.44575341) q[1];
x q[2];
rz(-0.85373803) q[3];
sx q[3];
rz(-1.8150363) q[3];
sx q[3];
rz(1.8239886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0539315) q[2];
sx q[2];
rz(-1.4940741) q[2];
sx q[2];
rz(1.9653448) q[2];
rz(-0.67115274) q[3];
sx q[3];
rz(-1.2495557) q[3];
sx q[3];
rz(-1.6254856) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072356) q[0];
sx q[0];
rz(-2.2800627) q[0];
sx q[0];
rz(-2.0980515) q[0];
rz(0.097298233) q[1];
sx q[1];
rz(-2.0847335) q[1];
sx q[1];
rz(2.0211438) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6098675) q[0];
sx q[0];
rz(-1.5104482) q[0];
sx q[0];
rz(-2.1670695) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78999577) q[2];
sx q[2];
rz(-2.1888791) q[2];
sx q[2];
rz(-2.2933985) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35353684) q[1];
sx q[1];
rz(-0.18419838) q[1];
sx q[1];
rz(0.26923979) q[1];
x q[2];
rz(-2.0529641) q[3];
sx q[3];
rz(-1.4034162) q[3];
sx q[3];
rz(-0.24576223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6843159) q[2];
sx q[2];
rz(-0.61351675) q[2];
sx q[2];
rz(1.3555917) q[2];
rz(-1.5654303) q[3];
sx q[3];
rz(-1.0836982) q[3];
sx q[3];
rz(2.1255597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23715699) q[0];
sx q[0];
rz(-2.5813527) q[0];
sx q[0];
rz(-0.78975633) q[0];
rz(2.4850028) q[1];
sx q[1];
rz(-2.0273392) q[1];
sx q[1];
rz(-0.56710342) q[1];
rz(0.48390735) q[2];
sx q[2];
rz(-0.96889062) q[2];
sx q[2];
rz(2.6775757) q[2];
rz(2.5593467) q[3];
sx q[3];
rz(-2.7979294) q[3];
sx q[3];
rz(2.9002849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
