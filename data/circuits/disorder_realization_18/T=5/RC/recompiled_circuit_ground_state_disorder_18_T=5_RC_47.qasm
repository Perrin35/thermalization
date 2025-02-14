OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9590149) q[0];
sx q[0];
rz(-2.9380517) q[0];
sx q[0];
rz(-1.8939053) q[0];
rz(1.3522476) q[1];
sx q[1];
rz(-0.69753733) q[1];
sx q[1];
rz(-0.69812671) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0188698) q[0];
sx q[0];
rz(-1.1636881) q[0];
sx q[0];
rz(-2.5886445) q[0];
x q[1];
rz(-0.45322541) q[2];
sx q[2];
rz(-1.5767323) q[2];
sx q[2];
rz(2.5763032) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7835603) q[1];
sx q[1];
rz(-0.81746819) q[1];
sx q[1];
rz(1.9897575) q[1];
rz(-0.014667347) q[3];
sx q[3];
rz(-1.8258182) q[3];
sx q[3];
rz(1.6153112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.84833604) q[2];
sx q[2];
rz(-1.9170599) q[2];
sx q[2];
rz(-0.87948925) q[2];
rz(0.73421156) q[3];
sx q[3];
rz(-2.0101533) q[3];
sx q[3];
rz(2.3122299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64049983) q[0];
sx q[0];
rz(-1.9556029) q[0];
sx q[0];
rz(-0.99047852) q[0];
rz(0.99539202) q[1];
sx q[1];
rz(-1.6973015) q[1];
sx q[1];
rz(-2.676414) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8841726) q[0];
sx q[0];
rz(-1.3883988) q[0];
sx q[0];
rz(0.42322961) q[0];
rz(-pi) q[1];
rz(2.1976333) q[2];
sx q[2];
rz(-1.8803681) q[2];
sx q[2];
rz(2.7322526) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.097660573) q[1];
sx q[1];
rz(-2.2273835) q[1];
sx q[1];
rz(0.93933479) q[1];
x q[2];
rz(1.7725132) q[3];
sx q[3];
rz(-1.3139825) q[3];
sx q[3];
rz(1.6981924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5140932) q[2];
sx q[2];
rz(-1.9072615) q[2];
sx q[2];
rz(2.9631183) q[2];
rz(0.56525362) q[3];
sx q[3];
rz(-1.3924799) q[3];
sx q[3];
rz(2.5843411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(2.8239215) q[0];
sx q[0];
rz(-1.3452106) q[0];
sx q[0];
rz(-2.4406261) q[0];
rz(2.7340381) q[1];
sx q[1];
rz(-1.2113672) q[1];
sx q[1];
rz(-2.0661381) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70301818) q[0];
sx q[0];
rz(-1.0319627) q[0];
sx q[0];
rz(1.9608905) q[0];
x q[1];
rz(-2.9161152) q[2];
sx q[2];
rz(-2.2144711) q[2];
sx q[2];
rz(0.89356092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9270529) q[1];
sx q[1];
rz(-0.42286122) q[1];
sx q[1];
rz(-1.8936669) q[1];
rz(-0.050450669) q[3];
sx q[3];
rz(-1.3063432) q[3];
sx q[3];
rz(1.3095625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.067523) q[2];
sx q[2];
rz(-1.3347551) q[2];
sx q[2];
rz(1.734181) q[2];
rz(-2.6635026) q[3];
sx q[3];
rz(-0.34308386) q[3];
sx q[3];
rz(-0.34496719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5673229) q[0];
sx q[0];
rz(-1.3936717) q[0];
sx q[0];
rz(1.0465013) q[0];
rz(-0.25431713) q[1];
sx q[1];
rz(-0.81552243) q[1];
sx q[1];
rz(0.3124803) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4506671) q[0];
sx q[0];
rz(-1.9798128) q[0];
sx q[0];
rz(3.0002007) q[0];
rz(-2.5529892) q[2];
sx q[2];
rz(-0.75102511) q[2];
sx q[2];
rz(-2.5680755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.65838012) q[1];
sx q[1];
rz(-1.8129895) q[1];
sx q[1];
rz(-0.32434742) q[1];
x q[2];
rz(-0.32799566) q[3];
sx q[3];
rz(-1.0803137) q[3];
sx q[3];
rz(-2.1759263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8326524) q[2];
sx q[2];
rz(-2.4231484) q[2];
sx q[2];
rz(-0.18860513) q[2];
rz(2.9077933) q[3];
sx q[3];
rz(-2.2457687) q[3];
sx q[3];
rz(-2.5578267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69498649) q[0];
sx q[0];
rz(-1.3124895) q[0];
sx q[0];
rz(-3.1332698) q[0];
rz(2.7024929) q[1];
sx q[1];
rz(-1.7726026) q[1];
sx q[1];
rz(-1.2459374) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6855658) q[0];
sx q[0];
rz(-1.7626806) q[0];
sx q[0];
rz(0.45996968) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14908423) q[2];
sx q[2];
rz(-1.8051762) q[2];
sx q[2];
rz(-2.5917376) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.627927) q[1];
sx q[1];
rz(-2.6509319) q[1];
sx q[1];
rz(-0.8497683) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28723081) q[3];
sx q[3];
rz(-0.93614791) q[3];
sx q[3];
rz(1.3040445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40450725) q[2];
sx q[2];
rz(-3.0948907) q[2];
sx q[2];
rz(-1.6339634) q[2];
rz(-2.5493933) q[3];
sx q[3];
rz(-1.3785572) q[3];
sx q[3];
rz(2.4018905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1208039) q[0];
sx q[0];
rz(-1.4494267) q[0];
sx q[0];
rz(2.8756496) q[0];
rz(1.2159011) q[1];
sx q[1];
rz(-0.78548702) q[1];
sx q[1];
rz(-1.9507834) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94534369) q[0];
sx q[0];
rz(-2.4085345) q[0];
sx q[0];
rz(-1.5668012) q[0];
rz(-pi) q[1];
rz(1.0180151) q[2];
sx q[2];
rz(-1.5399122) q[2];
sx q[2];
rz(-2.4190703) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.3336181) q[1];
sx q[1];
rz(-0.75389538) q[1];
sx q[1];
rz(3.1303166) q[1];
rz(-0.27482185) q[3];
sx q[3];
rz(-1.6259818) q[3];
sx q[3];
rz(1.7589057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33209458) q[2];
sx q[2];
rz(-1.7569434) q[2];
sx q[2];
rz(0.14796999) q[2];
rz(-2.68908) q[3];
sx q[3];
rz(-2.0544402) q[3];
sx q[3];
rz(-0.40542671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0026523503) q[0];
sx q[0];
rz(-1.2661221) q[0];
sx q[0];
rz(1.3462322) q[0];
rz(3.1345308) q[1];
sx q[1];
rz(-2.4738753) q[1];
sx q[1];
rz(-0.5074358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0327099) q[0];
sx q[0];
rz(-1.0531972) q[0];
sx q[0];
rz(0.93982307) q[0];
rz(1.9155923) q[2];
sx q[2];
rz(-1.7168247) q[2];
sx q[2];
rz(-0.70521077) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7120613) q[1];
sx q[1];
rz(-1.1768394) q[1];
sx q[1];
rz(-1.560657) q[1];
rz(1.8431013) q[3];
sx q[3];
rz(-1.57988) q[3];
sx q[3];
rz(0.42382012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.70600447) q[2];
sx q[2];
rz(-2.0861349) q[2];
sx q[2];
rz(0.17534193) q[2];
rz(2.7591779) q[3];
sx q[3];
rz(-1.9491111) q[3];
sx q[3];
rz(-2.6800938) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14313702) q[0];
sx q[0];
rz(-2.1303506) q[0];
sx q[0];
rz(1.3125516) q[0];
rz(3.0328499) q[1];
sx q[1];
rz(-0.60832945) q[1];
sx q[1];
rz(-2.463602) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16081339) q[0];
sx q[0];
rz(-0.49404544) q[0];
sx q[0];
rz(-1.52397) q[0];
x q[1];
rz(-0.021804734) q[2];
sx q[2];
rz(-0.74912375) q[2];
sx q[2];
rz(-2.5072012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7398518) q[1];
sx q[1];
rz(-1.051552) q[1];
sx q[1];
rz(0.3335764) q[1];
rz(-pi) q[2];
rz(2.1021712) q[3];
sx q[3];
rz(-2.2153004) q[3];
sx q[3];
rz(2.927305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2737736) q[2];
sx q[2];
rz(-1.5673693) q[2];
sx q[2];
rz(-1.1567953) q[2];
rz(-0.50446883) q[3];
sx q[3];
rz(-1.0320832) q[3];
sx q[3];
rz(-2.5696136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.2061763) q[0];
sx q[0];
rz(-1.6125866) q[0];
sx q[0];
rz(2.5231498) q[0];
rz(-2.2930324) q[1];
sx q[1];
rz(-2.6237374) q[1];
sx q[1];
rz(2.3458164) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2413809) q[0];
sx q[0];
rz(-2.6707895) q[0];
sx q[0];
rz(-1.6241252) q[0];
x q[1];
rz(0.82473849) q[2];
sx q[2];
rz(-1.5360263) q[2];
sx q[2];
rz(-0.29037133) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.56346169) q[1];
sx q[1];
rz(-1.8860779) q[1];
sx q[1];
rz(2.8540846) q[1];
x q[2];
rz(0.14988092) q[3];
sx q[3];
rz(-1.1544041) q[3];
sx q[3];
rz(0.021377953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.985118) q[2];
sx q[2];
rz(-1.5072344) q[2];
sx q[2];
rz(-1.9462684) q[2];
rz(-2.7965503) q[3];
sx q[3];
rz(-2.2796567) q[3];
sx q[3];
rz(0.8663469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
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
rz(1.7095551) q[0];
sx q[0];
rz(-2.9032752) q[0];
sx q[0];
rz(-3.0639783) q[0];
rz(1.2331102) q[1];
sx q[1];
rz(-1.0401007) q[1];
sx q[1];
rz(2.3495823) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.663765) q[0];
sx q[0];
rz(-0.91807967) q[0];
sx q[0];
rz(2.2774612) q[0];
rz(0.15969212) q[2];
sx q[2];
rz(-1.4174685) q[2];
sx q[2];
rz(-3.0042269) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35490552) q[1];
sx q[1];
rz(-0.49017683) q[1];
sx q[1];
rz(1.9016983) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11802063) q[3];
sx q[3];
rz(-0.80720085) q[3];
sx q[3];
rz(1.0420052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21976694) q[2];
sx q[2];
rz(-2.2302088) q[2];
sx q[2];
rz(2.0683973) q[2];
rz(-2.8940767) q[3];
sx q[3];
rz(-1.211834) q[3];
sx q[3];
rz(3.1246576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3794004) q[0];
sx q[0];
rz(-1.8397377) q[0];
sx q[0];
rz(2.4158438) q[0];
rz(0.87860592) q[1];
sx q[1];
rz(-2.286372) q[1];
sx q[1];
rz(1.1927037) q[1];
rz(2.0547619) q[2];
sx q[2];
rz(-1.3916343) q[2];
sx q[2];
rz(1.8050031) q[2];
rz(-1.5069275) q[3];
sx q[3];
rz(-0.72790925) q[3];
sx q[3];
rz(2.6016605) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
