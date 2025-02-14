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
rz(1.4662161) q[0];
sx q[0];
rz(5.3386547) q[0];
sx q[0];
rz(9.6261779) q[0];
rz(4.2138777) q[1];
sx q[1];
rz(3.9044851) q[1];
sx q[1];
rz(4.8621096) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9147707) q[0];
sx q[0];
rz(-0.96089191) q[0];
sx q[0];
rz(3.0746835) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6236963) q[2];
sx q[2];
rz(-0.96518789) q[2];
sx q[2];
rz(2.6014858) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0951251) q[1];
sx q[1];
rz(-0.78739843) q[1];
sx q[1];
rz(0.56708972) q[1];
rz(-pi) q[2];
rz(1.1789315) q[3];
sx q[3];
rz(-1.2614935) q[3];
sx q[3];
rz(-0.40180692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1066771) q[2];
sx q[2];
rz(-0.80433977) q[2];
sx q[2];
rz(2.8625281) q[2];
rz(-1.2532824) q[3];
sx q[3];
rz(-0.24975714) q[3];
sx q[3];
rz(-0.15160027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5308373) q[0];
sx q[0];
rz(-0.84722561) q[0];
sx q[0];
rz(-2.8952059) q[0];
rz(1.9465744) q[1];
sx q[1];
rz(-2.2476826) q[1];
sx q[1];
rz(-1.5066719) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6794066) q[0];
sx q[0];
rz(-2.6194318) q[0];
sx q[0];
rz(2.0798111) q[0];
rz(-0.042889281) q[2];
sx q[2];
rz(-1.873462) q[2];
sx q[2];
rz(-1.8406701) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40827258) q[1];
sx q[1];
rz(-0.36767861) q[1];
sx q[1];
rz(-1.7321083) q[1];
rz(-0.89119567) q[3];
sx q[3];
rz(-2.5053484) q[3];
sx q[3];
rz(-1.6980069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7404777) q[2];
sx q[2];
rz(-0.97492188) q[2];
sx q[2];
rz(-1.231989) q[2];
rz(1.1022107) q[3];
sx q[3];
rz(-3.1372742) q[3];
sx q[3];
rz(1.3365041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44322893) q[0];
sx q[0];
rz(-0.6969499) q[0];
sx q[0];
rz(2.2337636) q[0];
rz(2.5746131) q[1];
sx q[1];
rz(-2.1588529) q[1];
sx q[1];
rz(0.10279113) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55926052) q[0];
sx q[0];
rz(-2.1301422) q[0];
sx q[0];
rz(1.8190967) q[0];
rz(0.081016531) q[2];
sx q[2];
rz(-1.6747432) q[2];
sx q[2];
rz(-2.8250717) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6425848) q[1];
sx q[1];
rz(-0.73621589) q[1];
sx q[1];
rz(-2.5166488) q[1];
x q[2];
rz(-2.0538883) q[3];
sx q[3];
rz(-0.60988855) q[3];
sx q[3];
rz(2.1263378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28321442) q[2];
sx q[2];
rz(-0.80867043) q[2];
sx q[2];
rz(1.7041448) q[2];
rz(-0.39250675) q[3];
sx q[3];
rz(-2.0065353) q[3];
sx q[3];
rz(2.8431622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0482386) q[0];
sx q[0];
rz(-1.7788576) q[0];
sx q[0];
rz(-1.1887953) q[0];
rz(2.8554754) q[1];
sx q[1];
rz(-0.61758271) q[1];
sx q[1];
rz(-1.6292705) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85507369) q[0];
sx q[0];
rz(-0.60337702) q[0];
sx q[0];
rz(1.0666749) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9061162) q[2];
sx q[2];
rz(-2.309707) q[2];
sx q[2];
rz(-0.12140935) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.69763598) q[1];
sx q[1];
rz(-0.62503615) q[1];
sx q[1];
rz(-1.8915063) q[1];
x q[2];
rz(1.2028473) q[3];
sx q[3];
rz(-0.5371437) q[3];
sx q[3];
rz(0.14432553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7901223) q[2];
sx q[2];
rz(-2.0581547) q[2];
sx q[2];
rz(-2.4646087) q[2];
rz(1.2466768) q[3];
sx q[3];
rz(-1.2723943) q[3];
sx q[3];
rz(-1.5164794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93382728) q[0];
sx q[0];
rz(-2.8849869) q[0];
sx q[0];
rz(-0.024913464) q[0];
rz(-2.086153) q[1];
sx q[1];
rz(-0.98506227) q[1];
sx q[1];
rz(-0.28344646) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058539778) q[0];
sx q[0];
rz(-1.3136567) q[0];
sx q[0];
rz(-1.8399946) q[0];
x q[1];
rz(-1.3749429) q[2];
sx q[2];
rz(-1.699866) q[2];
sx q[2];
rz(-1.376525) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39549669) q[1];
sx q[1];
rz(-1.0701792) q[1];
sx q[1];
rz(-3.1067645) q[1];
rz(-pi) q[2];
rz(-0.33971649) q[3];
sx q[3];
rz(-0.92631683) q[3];
sx q[3];
rz(-3.0655053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.31263605) q[2];
sx q[2];
rz(-2.7241311) q[2];
sx q[2];
rz(2.8157595) q[2];
rz(1.2827986) q[3];
sx q[3];
rz(-1.157434) q[3];
sx q[3];
rz(-1.5939943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40389898) q[0];
sx q[0];
rz(-1.6588545) q[0];
sx q[0];
rz(0.37495908) q[0];
rz(-0.47863475) q[1];
sx q[1];
rz(-0.67239434) q[1];
sx q[1];
rz(0.37014827) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86287303) q[0];
sx q[0];
rz(-1.5476883) q[0];
sx q[0];
rz(1.5591168) q[0];
rz(-pi) q[1];
rz(-0.13932087) q[2];
sx q[2];
rz(-1.754154) q[2];
sx q[2];
rz(-0.56688165) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3016175) q[1];
sx q[1];
rz(-2.2711922) q[1];
sx q[1];
rz(-2.5216334) q[1];
rz(0.010995098) q[3];
sx q[3];
rz(-0.42783005) q[3];
sx q[3];
rz(1.7181991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39773539) q[2];
sx q[2];
rz(-0.19146679) q[2];
sx q[2];
rz(1.3810623) q[2];
rz(-1.0487652) q[3];
sx q[3];
rz(-1.7966813) q[3];
sx q[3];
rz(-0.63961187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25310707) q[0];
sx q[0];
rz(-0.83919224) q[0];
sx q[0];
rz(-2.2174368) q[0];
rz(2.2249075) q[1];
sx q[1];
rz(-1.0662096) q[1];
sx q[1];
rz(1.4013269) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8290212) q[0];
sx q[0];
rz(-1.628546) q[0];
sx q[0];
rz(2.6472378) q[0];
rz(-pi) q[1];
rz(2.9193814) q[2];
sx q[2];
rz(-0.66721081) q[2];
sx q[2];
rz(1.3363584) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3360916) q[1];
sx q[1];
rz(-0.84118836) q[1];
sx q[1];
rz(-1.7859687) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0592693) q[3];
sx q[3];
rz(-1.7632177) q[3];
sx q[3];
rz(-0.50979641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2283198) q[2];
sx q[2];
rz(-1.4618123) q[2];
sx q[2];
rz(1.224996) q[2];
rz(-3.1332968) q[3];
sx q[3];
rz(-3.1305997) q[3];
sx q[3];
rz(2.2744961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55815721) q[0];
sx q[0];
rz(-0.83034101) q[0];
sx q[0];
rz(-0.48026568) q[0];
rz(-0.14446124) q[1];
sx q[1];
rz(-2.2339349) q[1];
sx q[1];
rz(-2.2772148) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9821765) q[0];
sx q[0];
rz(-1.6173395) q[0];
sx q[0];
rz(-1.7345364) q[0];
rz(-pi) q[1];
rz(3.046918) q[2];
sx q[2];
rz(-1.5370768) q[2];
sx q[2];
rz(2.0168436) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4982953) q[1];
sx q[1];
rz(-1.3302696) q[1];
sx q[1];
rz(1.0161736) q[1];
rz(-pi) q[2];
rz(-1.8501353) q[3];
sx q[3];
rz(-1.8350826) q[3];
sx q[3];
rz(0.13388982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9353443) q[2];
sx q[2];
rz(-1.0103005) q[2];
sx q[2];
rz(-2.5065191) q[2];
rz(-2.5687929) q[3];
sx q[3];
rz(-1.7856995) q[3];
sx q[3];
rz(-0.87535453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11005814) q[0];
sx q[0];
rz(-0.77339554) q[0];
sx q[0];
rz(3.0173259) q[0];
rz(-0.36525137) q[1];
sx q[1];
rz(-1.4271913) q[1];
sx q[1];
rz(-1.7100547) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8339523) q[0];
sx q[0];
rz(-0.26276428) q[0];
sx q[0];
rz(0.74547474) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4489787) q[2];
sx q[2];
rz(-0.78348762) q[2];
sx q[2];
rz(0.88800511) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.86399779) q[1];
sx q[1];
rz(-0.91121948) q[1];
sx q[1];
rz(-1.9602174) q[1];
rz(-pi) q[2];
rz(-0.67691524) q[3];
sx q[3];
rz(-2.4999818) q[3];
sx q[3];
rz(0.42719242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26620904) q[2];
sx q[2];
rz(-0.92140809) q[2];
sx q[2];
rz(1.1727775) q[2];
rz(-1.734599) q[3];
sx q[3];
rz(-1.3318136) q[3];
sx q[3];
rz(1.2146568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5491972) q[0];
sx q[0];
rz(-2.0068491) q[0];
sx q[0];
rz(-0.59463516) q[0];
rz(1.0058962) q[1];
sx q[1];
rz(-1.9391831) q[1];
sx q[1];
rz(-0.88919052) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.517006) q[0];
sx q[0];
rz(-1.3505624) q[0];
sx q[0];
rz(0.051924719) q[0];
x q[1];
rz(-1.7125193) q[2];
sx q[2];
rz(-2.6772039) q[2];
sx q[2];
rz(1.022867) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8915705) q[1];
sx q[1];
rz(-1.6305939) q[1];
sx q[1];
rz(-1.476261) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5327644) q[3];
sx q[3];
rz(-1.0547148) q[3];
sx q[3];
rz(-0.79697679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1903926) q[2];
sx q[2];
rz(-1.9659646) q[2];
sx q[2];
rz(-0.54538837) q[2];
rz(-2.178318) q[3];
sx q[3];
rz(-1.1759718) q[3];
sx q[3];
rz(-1.4639328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6496898) q[0];
sx q[0];
rz(-0.90072537) q[0];
sx q[0];
rz(-1.8186722) q[0];
rz(2.2979965) q[1];
sx q[1];
rz(-1.0782764) q[1];
sx q[1];
rz(-1.2596399) q[1];
rz(-1.1314992) q[2];
sx q[2];
rz(-1.7514624) q[2];
sx q[2];
rz(-3.024586) q[2];
rz(2.4424408) q[3];
sx q[3];
rz(-0.71157645) q[3];
sx q[3];
rz(-0.012307766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
