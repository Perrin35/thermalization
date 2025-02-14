OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.71819031) q[0];
sx q[0];
rz(-3.1388404) q[0];
sx q[0];
rz(1.0116853) q[0];
rz(0.45971316) q[1];
sx q[1];
rz(4.4432321) q[1];
sx q[1];
rz(9.2530773) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9515686) q[0];
sx q[0];
rz(-2.5309238) q[0];
sx q[0];
rz(-0.14279731) q[0];
rz(0.5131716) q[2];
sx q[2];
rz(-2.9193239) q[2];
sx q[2];
rz(-1.5602416) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.57331177) q[1];
sx q[1];
rz(-2.3900095) q[1];
sx q[1];
rz(-2.4922396) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0875077) q[3];
sx q[3];
rz(-1.4776609) q[3];
sx q[3];
rz(-2.2786811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1823938) q[2];
sx q[2];
rz(-2.6540519) q[2];
sx q[2];
rz(1.4200776) q[2];
rz(-1.9567664) q[3];
sx q[3];
rz(-1.607837) q[3];
sx q[3];
rz(2.0737295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083953388) q[0];
sx q[0];
rz(-1.5204484) q[0];
sx q[0];
rz(2.6481096) q[0];
rz(2.3656942) q[1];
sx q[1];
rz(-2.6319365) q[1];
sx q[1];
rz(0.50481558) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017708721) q[0];
sx q[0];
rz(-1.0320745) q[0];
sx q[0];
rz(-0.81102844) q[0];
x q[1];
rz(0.52689441) q[2];
sx q[2];
rz(-1.829648) q[2];
sx q[2];
rz(0.35824725) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3146554) q[1];
sx q[1];
rz(-0.75410226) q[1];
sx q[1];
rz(-2.5376863) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8756469) q[3];
sx q[3];
rz(-1.6308745) q[3];
sx q[3];
rz(1.657965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33727553) q[2];
sx q[2];
rz(-1.3108459) q[2];
sx q[2];
rz(-0.47067019) q[2];
rz(0.63052952) q[3];
sx q[3];
rz(-1.0258521) q[3];
sx q[3];
rz(-1.4001018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077496342) q[0];
sx q[0];
rz(-3.016576) q[0];
sx q[0];
rz(0.65573829) q[0];
rz(-0.31133044) q[1];
sx q[1];
rz(-2.2475524) q[1];
sx q[1];
rz(-3.0701367) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0229683) q[0];
sx q[0];
rz(-0.35425348) q[0];
sx q[0];
rz(-1.7153165) q[0];
rz(-pi) q[1];
rz(-2.803903) q[2];
sx q[2];
rz(-1.9200385) q[2];
sx q[2];
rz(1.2472635) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.25117043) q[1];
sx q[1];
rz(-2.648472) q[1];
sx q[1];
rz(-1.7095196) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7903665) q[3];
sx q[3];
rz(-0.5595419) q[3];
sx q[3];
rz(2.9282687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91325703) q[2];
sx q[2];
rz(-1.5847289) q[2];
sx q[2];
rz(-3.081591) q[2];
rz(1.9289198) q[3];
sx q[3];
rz(-2.4433177) q[3];
sx q[3];
rz(-0.76446271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5948828) q[0];
sx q[0];
rz(-1.3285652) q[0];
sx q[0];
rz(2.5643964) q[0];
rz(2.3120841) q[1];
sx q[1];
rz(-2.0894158) q[1];
sx q[1];
rz(-0.83782354) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55623193) q[0];
sx q[0];
rz(-0.62667003) q[0];
sx q[0];
rz(-3.0100432) q[0];
x q[1];
rz(0.98060645) q[2];
sx q[2];
rz(-1.5204805) q[2];
sx q[2];
rz(-2.3286164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.29463331) q[1];
sx q[1];
rz(-0.75645743) q[1];
sx q[1];
rz(-1.0119757) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89452215) q[3];
sx q[3];
rz(-1.9274047) q[3];
sx q[3];
rz(1.5490393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.025658) q[2];
sx q[2];
rz(-1.6178774) q[2];
sx q[2];
rz(-1.9177829) q[2];
rz(2.6754248) q[3];
sx q[3];
rz(-0.40188447) q[3];
sx q[3];
rz(2.536072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25662988) q[0];
sx q[0];
rz(-1.4361359) q[0];
sx q[0];
rz(3.0404941) q[0];
rz(-2.7283607) q[1];
sx q[1];
rz(-2.7006472) q[1];
sx q[1];
rz(-2.0535927) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1285217) q[0];
sx q[0];
rz(-2.5422342) q[0];
sx q[0];
rz(-0.96512633) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1674983) q[2];
sx q[2];
rz(-2.3638862) q[2];
sx q[2];
rz(2.4571927) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1697547) q[1];
sx q[1];
rz(-2.7248451) q[1];
sx q[1];
rz(-0.93081148) q[1];
x q[2];
rz(-1.5002046) q[3];
sx q[3];
rz(-0.90259457) q[3];
sx q[3];
rz(-1.6833504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2130412) q[2];
sx q[2];
rz(-2.2553359) q[2];
sx q[2];
rz(1.586033) q[2];
rz(-2.0689615) q[3];
sx q[3];
rz(-1.5242679) q[3];
sx q[3];
rz(3.0090289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9685386) q[0];
sx q[0];
rz(-0.88468495) q[0];
sx q[0];
rz(1.435085) q[0];
rz(-0.75812078) q[1];
sx q[1];
rz(-2.4393612) q[1];
sx q[1];
rz(0.10163669) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91824965) q[0];
sx q[0];
rz(-2.0416284) q[0];
sx q[0];
rz(2.4001394) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36113895) q[2];
sx q[2];
rz(-1.1419347) q[2];
sx q[2];
rz(2.7782235) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1730301) q[1];
sx q[1];
rz(-0.30245879) q[1];
sx q[1];
rz(-2.7773501) q[1];
rz(-pi) q[2];
rz(2.6061344) q[3];
sx q[3];
rz(-2.1369918) q[3];
sx q[3];
rz(-0.2406075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8018084) q[2];
sx q[2];
rz(-1.1016176) q[2];
sx q[2];
rz(1.3327117) q[2];
rz(2.0594275) q[3];
sx q[3];
rz(-2.3887631) q[3];
sx q[3];
rz(-1.7180721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69671714) q[0];
sx q[0];
rz(-1.8901905) q[0];
sx q[0];
rz(-0.74309293) q[0];
rz(-0.7630868) q[1];
sx q[1];
rz(-2.5021195) q[1];
sx q[1];
rz(-2.2230164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7189499) q[0];
sx q[0];
rz(-0.026233679) q[0];
sx q[0];
rz(-1.8842741) q[0];
x q[1];
rz(0.3438216) q[2];
sx q[2];
rz(-2.0190416) q[2];
sx q[2];
rz(-0.78735414) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1700701) q[1];
sx q[1];
rz(-0.60107175) q[1];
sx q[1];
rz(2.5843589) q[1];
rz(2.4174446) q[3];
sx q[3];
rz(-1.482943) q[3];
sx q[3];
rz(-1.7334565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4172198) q[2];
sx q[2];
rz(-1.959274) q[2];
sx q[2];
rz(2.7057538) q[2];
rz(-3.0271652) q[3];
sx q[3];
rz(-0.2468214) q[3];
sx q[3];
rz(-2.54336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33335394) q[0];
sx q[0];
rz(-0.55753189) q[0];
sx q[0];
rz(-2.5575141) q[0];
rz(0.43560478) q[1];
sx q[1];
rz(-1.6807115) q[1];
sx q[1];
rz(-1.6216507) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0348957) q[0];
sx q[0];
rz(-0.55904065) q[0];
sx q[0];
rz(-0.99026545) q[0];
rz(-2.9233515) q[2];
sx q[2];
rz(-1.9710961) q[2];
sx q[2];
rz(-1.2117653) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4055109) q[1];
sx q[1];
rz(-0.20624948) q[1];
sx q[1];
rz(-3.095605) q[1];
x q[2];
rz(-2.5831843) q[3];
sx q[3];
rz(-2.8388151) q[3];
sx q[3];
rz(1.6807792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.074177563) q[2];
sx q[2];
rz(-1.7895074) q[2];
sx q[2];
rz(2.0254859) q[2];
rz(1.3793147) q[3];
sx q[3];
rz(-1.1032871) q[3];
sx q[3];
rz(3.0232159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24935687) q[0];
sx q[0];
rz(-0.07769575) q[0];
sx q[0];
rz(2.351601) q[0];
rz(-3.0780011) q[1];
sx q[1];
rz(-0.60706943) q[1];
sx q[1];
rz(-2.7008609) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1871619) q[0];
sx q[0];
rz(-1.8858028) q[0];
sx q[0];
rz(-0.31401547) q[0];
x q[1];
rz(-1.1665383) q[2];
sx q[2];
rz(-1.1102997) q[2];
sx q[2];
rz(1.3590036) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72183441) q[1];
sx q[1];
rz(-1.4065308) q[1];
sx q[1];
rz(1.0895776) q[1];
x q[2];
rz(3.1001631) q[3];
sx q[3];
rz(-1.3815615) q[3];
sx q[3];
rz(-1.3845058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.26238394) q[2];
sx q[2];
rz(-0.79664207) q[2];
sx q[2];
rz(-2.4764496) q[2];
rz(2.3696259) q[3];
sx q[3];
rz(-0.77163458) q[3];
sx q[3];
rz(1.6316679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7803698) q[0];
sx q[0];
rz(-2.4936115) q[0];
sx q[0];
rz(-1.004647) q[0];
rz(-1.7338344) q[1];
sx q[1];
rz(-1.6561457) q[1];
sx q[1];
rz(-1.2900603) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72899517) q[0];
sx q[0];
rz(-2.3069803) q[0];
sx q[0];
rz(0.33552977) q[0];
x q[1];
rz(1.2783706) q[2];
sx q[2];
rz(-2.5212581) q[2];
sx q[2];
rz(0.32493704) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.066034868) q[1];
sx q[1];
rz(-1.7820638) q[1];
sx q[1];
rz(-0.28055625) q[1];
x q[2];
rz(-2.9674888) q[3];
sx q[3];
rz(-0.977036) q[3];
sx q[3];
rz(-1.4623778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0066234) q[2];
sx q[2];
rz(-1.1639736) q[2];
sx q[2];
rz(2.8912365) q[2];
rz(-0.97744673) q[3];
sx q[3];
rz(-1.9340632) q[3];
sx q[3];
rz(1.6710963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.95494315) q[0];
sx q[0];
rz(-1.7087806) q[0];
sx q[0];
rz(1.9720672) q[0];
rz(-0.74465887) q[1];
sx q[1];
rz(-2.3458993) q[1];
sx q[1];
rz(0.69534272) q[1];
rz(2.9982243) q[2];
sx q[2];
rz(-1.5833387) q[2];
sx q[2];
rz(0.34064731) q[2];
rz(2.70784) q[3];
sx q[3];
rz(-2.0278145) q[3];
sx q[3];
rz(-0.47587004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
