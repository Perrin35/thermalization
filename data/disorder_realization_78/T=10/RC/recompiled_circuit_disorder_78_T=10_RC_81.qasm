OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5821563) q[0];
sx q[0];
rz(-0.59098935) q[0];
sx q[0];
rz(0.58340573) q[0];
rz(-0.18435873) q[1];
sx q[1];
rz(-2.1579722) q[1];
sx q[1];
rz(0.89259994) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8569782) q[0];
sx q[0];
rz(-1.0191139) q[0];
sx q[0];
rz(-0.3509699) q[0];
rz(-pi) q[1];
rz(1.0640261) q[2];
sx q[2];
rz(-0.96553409) q[2];
sx q[2];
rz(1.93047) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8436369) q[1];
sx q[1];
rz(-1.8814109) q[1];
sx q[1];
rz(2.288504) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6139908) q[3];
sx q[3];
rz(-2.3213263) q[3];
sx q[3];
rz(-0.84212069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9154174) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(1.9809451) q[2];
rz(2.9246269) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(2.0863566) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.083667) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(2.5657186) q[0];
rz(-1.2469762) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(-1.974568) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.089433) q[0];
sx q[0];
rz(-0.26425996) q[0];
sx q[0];
rz(1.7961851) q[0];
rz(-pi) q[1];
rz(0.24658792) q[2];
sx q[2];
rz(-1.0435259) q[2];
sx q[2];
rz(-1.2621244) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0547304) q[1];
sx q[1];
rz(-2.3332101) q[1];
sx q[1];
rz(3.0562835) q[1];
rz(2.4715273) q[3];
sx q[3];
rz(-1.9668005) q[3];
sx q[3];
rz(1.083064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.033826753) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(-0.21437422) q[2];
rz(3.068148) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(0.28081056) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6630994) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(1.2269155) q[0];
rz(-0.40027174) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(2.1267557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44714156) q[0];
sx q[0];
rz(-1.5146291) q[0];
sx q[0];
rz(-2.4980314) q[0];
rz(-pi) q[1];
rz(0.24200183) q[2];
sx q[2];
rz(-1.1761464) q[2];
sx q[2];
rz(1.5040656) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6940569) q[1];
sx q[1];
rz(-2.5390365) q[1];
sx q[1];
rz(1.7130997) q[1];
x q[2];
rz(-1.8869927) q[3];
sx q[3];
rz(-2.4305775) q[3];
sx q[3];
rz(2.7544114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0959452) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(0.8992368) q[2];
rz(-0.69747654) q[3];
sx q[3];
rz(-1.8557502) q[3];
sx q[3];
rz(-2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65524453) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(-0.43193257) q[0];
rz(-0.63255429) q[1];
sx q[1];
rz(-0.41699854) q[1];
sx q[1];
rz(2.5057709) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5084002) q[0];
sx q[0];
rz(-1.0754555) q[0];
sx q[0];
rz(0.20423996) q[0];
rz(1.7556778) q[2];
sx q[2];
rz(-0.73542483) q[2];
sx q[2];
rz(2.3787969) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1078474) q[1];
sx q[1];
rz(-0.81106942) q[1];
sx q[1];
rz(-1.8268405) q[1];
rz(-2.5084247) q[3];
sx q[3];
rz(-2.79106) q[3];
sx q[3];
rz(-1.4299973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2137961) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(2.4528743) q[2];
rz(-0.33411807) q[3];
sx q[3];
rz(-1.9709316) q[3];
sx q[3];
rz(-2.9978602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14389811) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(2.0671663) q[0];
rz(2.396446) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(2.863046) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6644088) q[0];
sx q[0];
rz(-1.9946949) q[0];
sx q[0];
rz(-0.82188481) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6804382) q[2];
sx q[2];
rz(-0.96193681) q[2];
sx q[2];
rz(-1.4102175) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.46036938) q[1];
sx q[1];
rz(-1.2037828) q[1];
sx q[1];
rz(-0.91194921) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3260818) q[3];
sx q[3];
rz(-1.1482571) q[3];
sx q[3];
rz(0.58327196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4429861) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(2.8473575) q[2];
rz(-3.0596628) q[3];
sx q[3];
rz(-2.6223845) q[3];
sx q[3];
rz(-0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(-0.0090573514) q[0];
rz(0.63502216) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(3.0335398) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8979643) q[0];
sx q[0];
rz(-1.122323) q[0];
sx q[0];
rz(-1.2178221) q[0];
rz(-pi) q[1];
rz(1.4322386) q[2];
sx q[2];
rz(-1.1922622) q[2];
sx q[2];
rz(-1.1059424) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3012645) q[1];
sx q[1];
rz(-2.5003308) q[1];
sx q[1];
rz(0.68323369) q[1];
rz(-pi) q[2];
x q[2];
rz(1.990854) q[3];
sx q[3];
rz(-0.58745158) q[3];
sx q[3];
rz(2.1849039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99284995) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(-1.2711058) q[2];
rz(3.0631915) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(-1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25061297) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(2.4839731) q[0];
rz(-1.3972067) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(-0.89362842) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2715627) q[0];
sx q[0];
rz(-1.8989519) q[0];
sx q[0];
rz(0.66332711) q[0];
x q[1];
rz(-1.0412752) q[2];
sx q[2];
rz(-0.5404226) q[2];
sx q[2];
rz(1.668001) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0221755) q[1];
sx q[1];
rz(-2.4319473) q[1];
sx q[1];
rz(-2.8303353) q[1];
rz(-pi) q[2];
rz(-1.4736389) q[3];
sx q[3];
rz(-2.328184) q[3];
sx q[3];
rz(1.1241476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39067337) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(-2.612109) q[2];
rz(2.6654065) q[3];
sx q[3];
rz(-1.3323077) q[3];
sx q[3];
rz(-0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3787518) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(1.09028) q[0];
rz(-0.11225637) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(1.9876678) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9268122) q[0];
sx q[0];
rz(-1.5367537) q[0];
sx q[0];
rz(-0.649931) q[0];
rz(-pi) q[1];
rz(1.4225142) q[2];
sx q[2];
rz(-2.1170756) q[2];
sx q[2];
rz(-0.15020457) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.47626074) q[1];
sx q[1];
rz(-1.6325103) q[1];
sx q[1];
rz(-2.3803821) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0941986) q[3];
sx q[3];
rz(-0.72924858) q[3];
sx q[3];
rz(0.13343982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(2.7611458) q[2];
rz(1.1278661) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(-1.5765223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34148759) q[0];
sx q[0];
rz(-1.2416168) q[0];
sx q[0];
rz(-2.5019116) q[0];
rz(-1.9027963) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(-1.170084) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6696475) q[0];
sx q[0];
rz(-1.8983316) q[0];
sx q[0];
rz(1.1963084) q[0];
rz(-pi) q[1];
rz(2.3652472) q[2];
sx q[2];
rz(-1.3753969) q[2];
sx q[2];
rz(-0.59567829) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.77800345) q[1];
sx q[1];
rz(-0.52007857) q[1];
sx q[1];
rz(2.7705454) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2421078) q[3];
sx q[3];
rz(-1.6791108) q[3];
sx q[3];
rz(1.9233821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8273948) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(-2.7588552) q[2];
rz(0.9283723) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-0.66463566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2255573) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(0.25892192) q[0];
rz(0.71031538) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(0.47992596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23586789) q[0];
sx q[0];
rz(-2.1972482) q[0];
sx q[0];
rz(0.41525526) q[0];
rz(-pi) q[1];
rz(-1.2725699) q[2];
sx q[2];
rz(-0.68694653) q[2];
sx q[2];
rz(-0.62703122) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5727947) q[1];
sx q[1];
rz(-0.63981445) q[1];
sx q[1];
rz(2.5261643) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94060937) q[3];
sx q[3];
rz(-2.3621231) q[3];
sx q[3];
rz(1.1704695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.84247983) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(1.3170362) q[2];
rz(-1.8995829) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(0.73808134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15923545) q[0];
sx q[0];
rz(-0.032082162) q[0];
sx q[0];
rz(-1.4627009) q[0];
rz(2.1622529) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(-0.32817763) q[2];
sx q[2];
rz(-2.6371418) q[2];
sx q[2];
rz(-0.20124659) q[2];
rz(0.76673037) q[3];
sx q[3];
rz(-2.7646716) q[3];
sx q[3];
rz(-0.9203831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];