OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7093631) q[0];
sx q[0];
rz(-2.1837283) q[0];
sx q[0];
rz(-0.14444484) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(-0.52539879) q[1];
sx q[1];
rz(-2.1638343) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9508096) q[0];
sx q[0];
rz(-1.3124183) q[0];
sx q[0];
rz(1.4825312) q[0];
x q[1];
rz(1.058504) q[2];
sx q[2];
rz(-2.5663178) q[2];
sx q[2];
rz(-3.1303867) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2797151) q[1];
sx q[1];
rz(-1.7126843) q[1];
sx q[1];
rz(0.14076294) q[1];
rz(-0.17840673) q[3];
sx q[3];
rz(-2.3462786) q[3];
sx q[3];
rz(0.91245302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4686761) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(0.93227512) q[2];
rz(-0.19876924) q[3];
sx q[3];
rz(-1.1105744) q[3];
sx q[3];
rz(0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4345877) q[0];
sx q[0];
rz(-2.2362237) q[0];
sx q[0];
rz(0.36112753) q[0];
rz(-1.4350285) q[1];
sx q[1];
rz(-1.3577434) q[1];
sx q[1];
rz(0.8180058) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084413962) q[0];
sx q[0];
rz(-1.955843) q[0];
sx q[0];
rz(-1.2731228) q[0];
rz(-pi) q[1];
rz(3.0188574) q[2];
sx q[2];
rz(-1.2018179) q[2];
sx q[2];
rz(0.69532794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.86707838) q[1];
sx q[1];
rz(-0.78650219) q[1];
sx q[1];
rz(2.2134476) q[1];
rz(-1.4217671) q[3];
sx q[3];
rz(-0.36747284) q[3];
sx q[3];
rz(3.0847103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8890185) q[2];
sx q[2];
rz(-2.694464) q[2];
sx q[2];
rz(-1.7209631) q[2];
rz(1.3160926) q[3];
sx q[3];
rz(-2.3836453) q[3];
sx q[3];
rz(-2.7533598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4017568) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(-0.87093583) q[0];
rz(0.31618205) q[1];
sx q[1];
rz(-2.8600287) q[1];
sx q[1];
rz(-0.2972163) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0270099) q[0];
sx q[0];
rz(-1.6342589) q[0];
sx q[0];
rz(0.40513904) q[0];
rz(-pi) q[1];
rz(0.88163968) q[2];
sx q[2];
rz(-2.214553) q[2];
sx q[2];
rz(1.8635441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9416562) q[1];
sx q[1];
rz(-1.7574649) q[1];
sx q[1];
rz(3.093064) q[1];
x q[2];
rz(2.475297) q[3];
sx q[3];
rz(-2.1994281) q[3];
sx q[3];
rz(-2.6578238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7685984) q[2];
sx q[2];
rz(-2.3186389) q[2];
sx q[2];
rz(-0.42759839) q[2];
rz(-1.1887431) q[3];
sx q[3];
rz(-2.5116428) q[3];
sx q[3];
rz(2.6141613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7052085) q[0];
sx q[0];
rz(-1.5773062) q[0];
sx q[0];
rz(-0.77392459) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(2.4598222) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1958774) q[0];
sx q[0];
rz(-1.3099075) q[0];
sx q[0];
rz(-2.9257141) q[0];
x q[1];
rz(2.9245124) q[2];
sx q[2];
rz(-2.1639369) q[2];
sx q[2];
rz(-0.3074239) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6097648) q[1];
sx q[1];
rz(-1.0351287) q[1];
sx q[1];
rz(0.51171724) q[1];
rz(-pi) q[2];
rz(-0.34337266) q[3];
sx q[3];
rz(-1.9150754) q[3];
sx q[3];
rz(-1.6791277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.008808) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(2.183389) q[2];
rz(-1.0664252) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(-0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5565857) q[0];
sx q[0];
rz(-1.7946694) q[0];
sx q[0];
rz(2.5812896) q[0];
rz(-0.99984461) q[1];
sx q[1];
rz(-2.9381349) q[1];
sx q[1];
rz(1.6220185) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1006267) q[0];
sx q[0];
rz(-2.0140411) q[0];
sx q[0];
rz(2.6020781) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72603307) q[2];
sx q[2];
rz(-0.75349977) q[2];
sx q[2];
rz(-1.7475278) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2677742) q[1];
sx q[1];
rz(-2.3582637) q[1];
sx q[1];
rz(2.8591213) q[1];
rz(0.95789692) q[3];
sx q[3];
rz(-1.4688204) q[3];
sx q[3];
rz(-1.0979872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4524298) q[2];
sx q[2];
rz(-2.5897557) q[2];
sx q[2];
rz(-0.2229283) q[2];
rz(0.034742268) q[3];
sx q[3];
rz(-1.3870753) q[3];
sx q[3];
rz(-3.070014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.3361622) q[0];
sx q[0];
rz(-2.8283089) q[0];
sx q[0];
rz(-2.0715332) q[0];
rz(1.7806212) q[1];
sx q[1];
rz(-2.762251) q[1];
sx q[1];
rz(1.8575352) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.703008) q[0];
sx q[0];
rz(-0.87235886) q[0];
sx q[0];
rz(-2.3914778) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.861704) q[2];
sx q[2];
rz(-1.1942689) q[2];
sx q[2];
rz(1.0371475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9575427) q[1];
sx q[1];
rz(-1.2831389) q[1];
sx q[1];
rz(2.6581453) q[1];
x q[2];
rz(1.0345801) q[3];
sx q[3];
rz(-1.3048733) q[3];
sx q[3];
rz(-1.769161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6665035) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(0.63759032) q[2];
rz(-2.3049138) q[3];
sx q[3];
rz(-2.0357318) q[3];
sx q[3];
rz(1.7975413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56629431) q[0];
sx q[0];
rz(-2.7116382) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(-2.7138846) q[1];
sx q[1];
rz(-1.5274915) q[1];
sx q[1];
rz(2.2033851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095031247) q[0];
sx q[0];
rz(-1.1633658) q[0];
sx q[0];
rz(2.808232) q[0];
rz(-pi) q[1];
rz(2.5222048) q[2];
sx q[2];
rz(-0.833138) q[2];
sx q[2];
rz(-2.7433861) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6066305) q[1];
sx q[1];
rz(-2.3772847) q[1];
sx q[1];
rz(2.0984089) q[1];
rz(-pi) q[2];
rz(-0.70721831) q[3];
sx q[3];
rz(-1.3361317) q[3];
sx q[3];
rz(2.9747687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2618835) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(-1.7555457) q[2];
rz(-1.8188247) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(0.02903207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26738527) q[0];
sx q[0];
rz(-2.8081942) q[0];
sx q[0];
rz(-1.7077131) q[0];
rz(1.2738312) q[1];
sx q[1];
rz(-2.0055983) q[1];
sx q[1];
rz(-0.83126718) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600663) q[0];
sx q[0];
rz(-0.43018451) q[0];
sx q[0];
rz(2.7823886) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7449042) q[2];
sx q[2];
rz(-2.499352) q[2];
sx q[2];
rz(-1.0345392) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9030994) q[1];
sx q[1];
rz(-2.0598754) q[1];
sx q[1];
rz(-2.3563983) q[1];
rz(0.76922272) q[3];
sx q[3];
rz(-2.8260494) q[3];
sx q[3];
rz(3.0358918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5780118) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(-0.9643628) q[2];
rz(1.9780805) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(0.65892974) q[3];
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
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7678541) q[0];
sx q[0];
rz(-2.4066194) q[0];
sx q[0];
rz(-2.1642165) q[0];
rz(1.7550229) q[1];
sx q[1];
rz(-1.8354548) q[1];
sx q[1];
rz(1.1057373) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2290105) q[0];
sx q[0];
rz(-1.6874896) q[0];
sx q[0];
rz(0.44793655) q[0];
rz(-2.2374723) q[2];
sx q[2];
rz(-1.6534001) q[2];
sx q[2];
rz(0.77043515) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8399664) q[1];
sx q[1];
rz(-1.3857538) q[1];
sx q[1];
rz(1.3355096) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2826142) q[3];
sx q[3];
rz(-1.2328706) q[3];
sx q[3];
rz(1.5087138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5332807) q[2];
sx q[2];
rz(-0.38485843) q[2];
sx q[2];
rz(2.9917955) q[2];
rz(-1.7685361) q[3];
sx q[3];
rz(-1.405973) q[3];
sx q[3];
rz(1.8201374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49366429) q[0];
sx q[0];
rz(-0.51769185) q[0];
sx q[0];
rz(0.1272442) q[0];
rz(1.6607025) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(0.1677992) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084899336) q[0];
sx q[0];
rz(-3.1057682) q[0];
sx q[0];
rz(-2.5762659) q[0];
x q[1];
rz(-2.5334353) q[2];
sx q[2];
rz(-1.4434575) q[2];
sx q[2];
rz(1.5437768) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64393109) q[1];
sx q[1];
rz(-1.5776458) q[1];
sx q[1];
rz(0.58845206) q[1];
x q[2];
rz(2.3778902) q[3];
sx q[3];
rz(-0.16994952) q[3];
sx q[3];
rz(0.22596879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.026713513) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(2.3804469) q[2];
rz(-0.090027697) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(-2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54820838) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(2.7535915) q[1];
sx q[1];
rz(-1.3996268) q[1];
sx q[1];
rz(-0.7849801) q[1];
rz(2.3315196) q[2];
sx q[2];
rz(-2.0740866) q[2];
sx q[2];
rz(-1.496051) q[2];
rz(1.2896982) q[3];
sx q[3];
rz(-0.55681183) q[3];
sx q[3];
rz(-1.1435215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
