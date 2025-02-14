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
rz(-1.7419185) q[0];
sx q[0];
rz(-1.4631441) q[0];
sx q[0];
rz(1.4611257) q[0];
rz(-0.48179102) q[1];
sx q[1];
rz(-0.81463373) q[1];
sx q[1];
rz(-2.2021267) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0793014) q[0];
sx q[0];
rz(-1.5271375) q[0];
sx q[0];
rz(0.0096304499) q[0];
rz(-pi) q[1];
rz(-1.9389802) q[2];
sx q[2];
rz(-2.0716616) q[2];
sx q[2];
rz(3.0291968) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0306326) q[1];
sx q[1];
rz(-0.041555066) q[1];
sx q[1];
rz(-1.3135733) q[1];
rz(1.788673) q[3];
sx q[3];
rz(-1.1718201) q[3];
sx q[3];
rz(0.50673317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7751094) q[2];
sx q[2];
rz(-1.6747687) q[2];
sx q[2];
rz(-0.62466204) q[2];
rz(2.1810253) q[3];
sx q[3];
rz(-1.4550236) q[3];
sx q[3];
rz(0.87489405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7333154) q[0];
sx q[0];
rz(-0.69077078) q[0];
sx q[0];
rz(-0.1917924) q[0];
rz(0.57418144) q[1];
sx q[1];
rz(-1.6185113) q[1];
sx q[1];
rz(-0.40723732) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2977266) q[0];
sx q[0];
rz(-2.9653984) q[0];
sx q[0];
rz(-1.5017444) q[0];
rz(-1.919433) q[2];
sx q[2];
rz(-0.18770175) q[2];
sx q[2];
rz(2.3003464) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.25000962) q[1];
sx q[1];
rz(-2.574157) q[1];
sx q[1];
rz(-2.7944348) q[1];
x q[2];
rz(1.7329747) q[3];
sx q[3];
rz(-0.50807488) q[3];
sx q[3];
rz(-1.8014419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9521956) q[2];
sx q[2];
rz(-1.962683) q[2];
sx q[2];
rz(-0.65323812) q[2];
rz(2.2630528) q[3];
sx q[3];
rz(-1.0983175) q[3];
sx q[3];
rz(3.0285335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1427796) q[0];
sx q[0];
rz(-1.7072562) q[0];
sx q[0];
rz(1.0502195) q[0];
rz(2.5255919) q[1];
sx q[1];
rz(-1.417825) q[1];
sx q[1];
rz(2.2284257) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28024689) q[0];
sx q[0];
rz(-1.4644091) q[0];
sx q[0];
rz(1.4004571) q[0];
rz(1.9550858) q[2];
sx q[2];
rz(-2.8336513) q[2];
sx q[2];
rz(-0.66288391) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.84308594) q[1];
sx q[1];
rz(-1.8131327) q[1];
sx q[1];
rz(-2.0451115) q[1];
rz(1.1230041) q[3];
sx q[3];
rz(-0.7402484) q[3];
sx q[3];
rz(-2.972176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22127613) q[2];
sx q[2];
rz(-2.2344799) q[2];
sx q[2];
rz(0.059180666) q[2];
rz(-2.3024998) q[3];
sx q[3];
rz(-1.2181166) q[3];
sx q[3];
rz(-2.2742417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7729823) q[0];
sx q[0];
rz(-2.1324069) q[0];
sx q[0];
rz(0.34977812) q[0];
rz(-1.7971669) q[1];
sx q[1];
rz(-2.5549922) q[1];
sx q[1];
rz(0.1127359) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7736334) q[0];
sx q[0];
rz(-1.2293574) q[0];
sx q[0];
rz(2.6608047) q[0];
rz(-pi) q[1];
rz(-2.2305528) q[2];
sx q[2];
rz(-1.8000812) q[2];
sx q[2];
rz(-2.2509172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1904691) q[1];
sx q[1];
rz(-0.72812041) q[1];
sx q[1];
rz(1.1090167) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68513667) q[3];
sx q[3];
rz(-2.1129012) q[3];
sx q[3];
rz(2.477081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.747644) q[2];
sx q[2];
rz(-1.5851721) q[2];
sx q[2];
rz(-0.057223884) q[2];
rz(-1.8852437) q[3];
sx q[3];
rz(-0.57259721) q[3];
sx q[3];
rz(-0.73067874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2644065) q[0];
sx q[0];
rz(-1.1671966) q[0];
sx q[0];
rz(-0.71006829) q[0];
rz(-0.26142985) q[1];
sx q[1];
rz(-1.556004) q[1];
sx q[1];
rz(0.95103055) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79890811) q[0];
sx q[0];
rz(-2.5747188) q[0];
sx q[0];
rz(2.5004412) q[0];
x q[1];
rz(2.0567953) q[2];
sx q[2];
rz(-0.99691454) q[2];
sx q[2];
rz(-1.4950372) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3465865) q[1];
sx q[1];
rz(-1.7730496) q[1];
sx q[1];
rz(0.70960416) q[1];
rz(-pi) q[2];
rz(-0.099154648) q[3];
sx q[3];
rz(-2.2014849) q[3];
sx q[3];
rz(-0.48336467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.49932536) q[2];
sx q[2];
rz(-0.32411164) q[2];
sx q[2];
rz(-0.30936852) q[2];
rz(1.0950836) q[3];
sx q[3];
rz(-1.1284004) q[3];
sx q[3];
rz(-0.40422082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9139022) q[0];
sx q[0];
rz(-2.9588283) q[0];
sx q[0];
rz(2.2392654) q[0];
rz(3.0790216) q[1];
sx q[1];
rz(-1.2689271) q[1];
sx q[1];
rz(1.8960457) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6645836) q[0];
sx q[0];
rz(-2.4467108) q[0];
sx q[0];
rz(-0.24212153) q[0];
x q[1];
rz(1.6643413) q[2];
sx q[2];
rz(-2.5999843) q[2];
sx q[2];
rz(-0.53047413) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.77038232) q[1];
sx q[1];
rz(-0.66714811) q[1];
sx q[1];
rz(-0.70828153) q[1];
x q[2];
rz(-2.5596722) q[3];
sx q[3];
rz(-0.76184638) q[3];
sx q[3];
rz(0.91530734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41783276) q[2];
sx q[2];
rz(-1.87169) q[2];
sx q[2];
rz(-1.343824) q[2];
rz(1.2835245) q[3];
sx q[3];
rz(-2.696974) q[3];
sx q[3];
rz(2.4166687) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.035324) q[0];
sx q[0];
rz(-3.1043053) q[0];
sx q[0];
rz(-2.7653747) q[0];
rz(0.47209921) q[1];
sx q[1];
rz(-2.4596228) q[1];
sx q[1];
rz(-0.45904407) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5556407) q[0];
sx q[0];
rz(-1.9919073) q[0];
sx q[0];
rz(0.79035203) q[0];
rz(-pi) q[1];
rz(-0.44957513) q[2];
sx q[2];
rz(-1.7284942) q[2];
sx q[2];
rz(1.1166935) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9487293) q[1];
sx q[1];
rz(-1.7336083) q[1];
sx q[1];
rz(1.5019187) q[1];
x q[2];
rz(2.543942) q[3];
sx q[3];
rz(-1.5603746) q[3];
sx q[3];
rz(-2.8437993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.072319357) q[2];
sx q[2];
rz(-0.71162144) q[2];
sx q[2];
rz(-2.7834328) q[2];
rz(-0.88851309) q[3];
sx q[3];
rz(-1.0212967) q[3];
sx q[3];
rz(0.93069881) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096980378) q[0];
sx q[0];
rz(-2.8428069) q[0];
sx q[0];
rz(2.1646747) q[0];
rz(2.3878429) q[1];
sx q[1];
rz(-1.6981373) q[1];
sx q[1];
rz(1.6389729) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.174265) q[0];
sx q[0];
rz(-2.2696724) q[0];
sx q[0];
rz(-1.4682659) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6303275) q[2];
sx q[2];
rz(-2.3663372) q[2];
sx q[2];
rz(1.1182736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.21587196) q[1];
sx q[1];
rz(-1.2684857) q[1];
sx q[1];
rz(-0.69958256) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8863932) q[3];
sx q[3];
rz(-0.64794174) q[3];
sx q[3];
rz(1.2256988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1961907) q[2];
sx q[2];
rz(-1.6795029) q[2];
sx q[2];
rz(-2.2958344) q[2];
rz(-0.78019199) q[3];
sx q[3];
rz(-1.7834241) q[3];
sx q[3];
rz(-0.63854027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3188401) q[0];
sx q[0];
rz(-1.0718811) q[0];
sx q[0];
rz(0.34291294) q[0];
rz(-2.137939) q[1];
sx q[1];
rz(-1.2316848) q[1];
sx q[1];
rz(-0.71663219) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9166853) q[0];
sx q[0];
rz(-1.9193847) q[0];
sx q[0];
rz(0.7734359) q[0];
rz(2.5942973) q[2];
sx q[2];
rz(-0.5623874) q[2];
sx q[2];
rz(-0.15531103) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1799419) q[1];
sx q[1];
rz(-1.838521) q[1];
sx q[1];
rz(2.5336392) q[1];
x q[2];
rz(1.9699205) q[3];
sx q[3];
rz(-0.93825996) q[3];
sx q[3];
rz(2.5312531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0922682) q[2];
sx q[2];
rz(-1.6246139) q[2];
sx q[2];
rz(2.1410904) q[2];
rz(-2.3246824) q[3];
sx q[3];
rz(-1.2714081) q[3];
sx q[3];
rz(2.7774155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3889121) q[0];
sx q[0];
rz(-2.2081544) q[0];
sx q[0];
rz(0.77891427) q[0];
rz(0.79175103) q[1];
sx q[1];
rz(-1.9154895) q[1];
sx q[1];
rz(2.7955999) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35157555) q[0];
sx q[0];
rz(-1.6540961) q[0];
sx q[0];
rz(0.44805713) q[0];
rz(-pi) q[1];
x q[1];
rz(0.012076186) q[2];
sx q[2];
rz(-1.2674971) q[2];
sx q[2];
rz(0.30840087) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78948632) q[1];
sx q[1];
rz(-0.54692422) q[1];
sx q[1];
rz(1.5032306) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9729887) q[3];
sx q[3];
rz(-2.3762083) q[3];
sx q[3];
rz(-0.51494288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.148823) q[2];
sx q[2];
rz(-1.6042446) q[2];
sx q[2];
rz(-2.0796622) q[2];
rz(0.81260931) q[3];
sx q[3];
rz(-1.142623) q[3];
sx q[3];
rz(-0.44212166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84520311) q[0];
sx q[0];
rz(-2.3136105) q[0];
sx q[0];
rz(-1.4985341) q[0];
rz(-2.8614112) q[1];
sx q[1];
rz(-1.9277086) q[1];
sx q[1];
rz(2.3452506) q[1];
rz(-0.40899271) q[2];
sx q[2];
rz(-1.6074642) q[2];
sx q[2];
rz(1.4516914) q[2];
rz(2.0186043) q[3];
sx q[3];
rz(-0.79389865) q[3];
sx q[3];
rz(0.47108847) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
