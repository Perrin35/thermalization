OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3341137) q[0];
sx q[0];
rz(-2.4456094) q[0];
sx q[0];
rz(0.88710436) q[0];
rz(-0.71169418) q[1];
sx q[1];
rz(-1.0625755) q[1];
sx q[1];
rz(-0.33023155) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6905211) q[0];
sx q[0];
rz(-2.1981648) q[0];
sx q[0];
rz(-0.55212195) q[0];
rz(-pi) q[1];
rz(-1.2057958) q[2];
sx q[2];
rz(-2.096039) q[2];
sx q[2];
rz(-2.989632) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5042233) q[1];
sx q[1];
rz(-1.0801225) q[1];
sx q[1];
rz(1.4549903) q[1];
x q[2];
rz(0.078465538) q[3];
sx q[3];
rz(-2.0002504) q[3];
sx q[3];
rz(2.0450908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3868788) q[2];
sx q[2];
rz(-0.42749307) q[2];
sx q[2];
rz(1.0270366) q[2];
rz(2.2586281) q[3];
sx q[3];
rz(-1.3436147) q[3];
sx q[3];
rz(0.55356717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3565732) q[0];
sx q[0];
rz(-3.0632601) q[0];
sx q[0];
rz(-2.8739492) q[0];
rz(-1.74125) q[1];
sx q[1];
rz(-2.5282271) q[1];
sx q[1];
rz(0.5563446) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69795495) q[0];
sx q[0];
rz(-0.55715269) q[0];
sx q[0];
rz(-2.055026) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.914996) q[2];
sx q[2];
rz(-2.2453893) q[2];
sx q[2];
rz(3.0795902) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.90168629) q[1];
sx q[1];
rz(-2.1237897) q[1];
sx q[1];
rz(-0.4497952) q[1];
rz(-pi) q[2];
rz(0.58014262) q[3];
sx q[3];
rz(-2.1506967) q[3];
sx q[3];
rz(2.4471791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9050425) q[2];
sx q[2];
rz(-1.4507797) q[2];
sx q[2];
rz(1.2518008) q[2];
rz(-1.2647789) q[3];
sx q[3];
rz(-2.4301961) q[3];
sx q[3];
rz(1.0714162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4816137) q[0];
sx q[0];
rz(-2.5857506) q[0];
sx q[0];
rz(-2.3023093) q[0];
rz(1.6370157) q[1];
sx q[1];
rz(-1.4312276) q[1];
sx q[1];
rz(1.379871) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8913075) q[0];
sx q[0];
rz(-1.6943185) q[0];
sx q[0];
rz(3.059292) q[0];
x q[1];
rz(0.21442757) q[2];
sx q[2];
rz(-1.0357719) q[2];
sx q[2];
rz(0.24033879) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7127925) q[1];
sx q[1];
rz(-1.2512491) q[1];
sx q[1];
rz(2.0682813) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1228509) q[3];
sx q[3];
rz(-2.0716487) q[3];
sx q[3];
rz(-0.93397442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.442349) q[2];
sx q[2];
rz(-1.9244497) q[2];
sx q[2];
rz(-0.71845636) q[2];
rz(0.62026223) q[3];
sx q[3];
rz(-2.2060427) q[3];
sx q[3];
rz(0.81234318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45850596) q[0];
sx q[0];
rz(-1.2309256) q[0];
sx q[0];
rz(-1.7339535) q[0];
rz(-0.55157026) q[1];
sx q[1];
rz(-0.11081129) q[1];
sx q[1];
rz(-1.3161906) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7421417) q[0];
sx q[0];
rz(-0.67614284) q[0];
sx q[0];
rz(-2.0318982) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1291601) q[2];
sx q[2];
rz(-1.2913449) q[2];
sx q[2];
rz(-0.81571992) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.516663) q[1];
sx q[1];
rz(-1.9716814) q[1];
sx q[1];
rz(-1.4275803) q[1];
rz(-pi) q[2];
rz(0.67263453) q[3];
sx q[3];
rz(-2.1676873) q[3];
sx q[3];
rz(-1.2850645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5482285) q[2];
sx q[2];
rz(-2.7061988) q[2];
sx q[2];
rz(-1.0080053) q[2];
rz(2.8625782) q[3];
sx q[3];
rz(-0.81441003) q[3];
sx q[3];
rz(-2.6830955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8654883) q[0];
sx q[0];
rz(-1.3300329) q[0];
sx q[0];
rz(2.8651067) q[0];
rz(0.28929389) q[1];
sx q[1];
rz(-1.9464867) q[1];
sx q[1];
rz(2.8791265) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97476995) q[0];
sx q[0];
rz(-1.747805) q[0];
sx q[0];
rz(2.9357852) q[0];
rz(-pi) q[1];
rz(3.1308181) q[2];
sx q[2];
rz(-0.97948217) q[2];
sx q[2];
rz(0.39515218) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5804606) q[1];
sx q[1];
rz(-1.9101686) q[1];
sx q[1];
rz(0.58402337) q[1];
rz(-pi) q[2];
rz(-0.57837242) q[3];
sx q[3];
rz(-1.4422073) q[3];
sx q[3];
rz(-1.2728387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.7485973) q[2];
sx q[2];
rz(-2.5454919) q[2];
sx q[2];
rz(-1.3470915) q[2];
rz(0.10284452) q[3];
sx q[3];
rz(-1.9072073) q[3];
sx q[3];
rz(1.5170521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0031925072) q[0];
sx q[0];
rz(-1.6380558) q[0];
sx q[0];
rz(0.060977161) q[0];
rz(-1.2334476) q[1];
sx q[1];
rz(-1.3950709) q[1];
sx q[1];
rz(-1.5256418) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73404445) q[0];
sx q[0];
rz(-1.5015232) q[0];
sx q[0];
rz(2.6114527) q[0];
rz(-pi) q[1];
rz(1.5360918) q[2];
sx q[2];
rz(-2.3922709) q[2];
sx q[2];
rz(2.9108436) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9445584) q[1];
sx q[1];
rz(-1.4602293) q[1];
sx q[1];
rz(-1.8227897) q[1];
rz(-pi) q[2];
rz(-0.16651972) q[3];
sx q[3];
rz(-0.92574471) q[3];
sx q[3];
rz(-1.695961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47152758) q[2];
sx q[2];
rz(-1.7933041) q[2];
sx q[2];
rz(1.8111551) q[2];
rz(2.0165675) q[3];
sx q[3];
rz(-0.87441134) q[3];
sx q[3];
rz(-0.037954656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9576981) q[0];
sx q[0];
rz(-1.4426458) q[0];
sx q[0];
rz(-2.884602) q[0];
rz(-0.43486241) q[1];
sx q[1];
rz(-2.3915274) q[1];
sx q[1];
rz(-0.90352568) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5241476) q[0];
sx q[0];
rz(-0.092369583) q[0];
sx q[0];
rz(-1.6156107) q[0];
rz(-1.7062188) q[2];
sx q[2];
rz(-2.3122182) q[2];
sx q[2];
rz(-0.00450762) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3147887) q[1];
sx q[1];
rz(-1.5233524) q[1];
sx q[1];
rz(-3.080009) q[1];
rz(-pi) q[2];
rz(0.89374505) q[3];
sx q[3];
rz(-0.66679685) q[3];
sx q[3];
rz(0.1689942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8187108) q[2];
sx q[2];
rz(-1.633753) q[2];
sx q[2];
rz(1.4212849) q[2];
rz(0.70982248) q[3];
sx q[3];
rz(-2.0028508) q[3];
sx q[3];
rz(0.53409725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98599559) q[0];
sx q[0];
rz(-2.7105712) q[0];
sx q[0];
rz(-2.0565597) q[0];
rz(2.494508) q[1];
sx q[1];
rz(-1.1754464) q[1];
sx q[1];
rz(-3.0771902) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2213106) q[0];
sx q[0];
rz(-0.94641173) q[0];
sx q[0];
rz(1.1657752) q[0];
x q[1];
rz(2.3008806) q[2];
sx q[2];
rz(-0.86879769) q[2];
sx q[2];
rz(0.15499767) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6013424) q[1];
sx q[1];
rz(-0.86886084) q[1];
sx q[1];
rz(2.1871845) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2256718) q[3];
sx q[3];
rz(-0.43817156) q[3];
sx q[3];
rz(-0.34153356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25962466) q[2];
sx q[2];
rz(-2.7590064) q[2];
sx q[2];
rz(0.043370334) q[2];
rz(2.4671593) q[3];
sx q[3];
rz(-1.3573656) q[3];
sx q[3];
rz(1.1446713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44716537) q[0];
sx q[0];
rz(-2.5612216) q[0];
sx q[0];
rz(0.32661435) q[0];
rz(-2.0475552) q[1];
sx q[1];
rz(-2.2973165) q[1];
sx q[1];
rz(2.6577267) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.992386) q[0];
sx q[0];
rz(-1.0629553) q[0];
sx q[0];
rz(-2.3293428) q[0];
rz(2.7174453) q[2];
sx q[2];
rz(-1.7778991) q[2];
sx q[2];
rz(-2.066156) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6623508) q[1];
sx q[1];
rz(-1.7852946) q[1];
sx q[1];
rz(-2.3420326) q[1];
x q[2];
rz(-0.30535474) q[3];
sx q[3];
rz(-0.66514665) q[3];
sx q[3];
rz(1.0843474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4491552) q[2];
sx q[2];
rz(-2.1026976) q[2];
sx q[2];
rz(3.0628487) q[2];
rz(0.76830831) q[3];
sx q[3];
rz(-1.2945622) q[3];
sx q[3];
rz(2.9505742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67097265) q[0];
sx q[0];
rz(-2.5363531) q[0];
sx q[0];
rz(-3.0249) q[0];
rz(-2.281588) q[1];
sx q[1];
rz(-1.8000894) q[1];
sx q[1];
rz(-2.2681627) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1310932) q[0];
sx q[0];
rz(-2.2981055) q[0];
sx q[0];
rz(-2.8167679) q[0];
x q[1];
rz(-1.248718) q[2];
sx q[2];
rz(-2.0587741) q[2];
sx q[2];
rz(-2.8805594) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22120093) q[1];
sx q[1];
rz(-1.6947692) q[1];
sx q[1];
rz(-2.6392379) q[1];
rz(2.5301873) q[3];
sx q[3];
rz(-1.5152001) q[3];
sx q[3];
rz(1.9724595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7200835) q[2];
sx q[2];
rz(-2.8953711) q[2];
sx q[2];
rz(-0.949238) q[2];
rz(1.2468437) q[3];
sx q[3];
rz(-1.6722164) q[3];
sx q[3];
rz(-0.59103549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.128189) q[0];
sx q[0];
rz(-2.1407776) q[0];
sx q[0];
rz(1.5110973) q[0];
rz(-0.36698256) q[1];
sx q[1];
rz(-1.783168) q[1];
sx q[1];
rz(2.5617243) q[1];
rz(0.30324528) q[2];
sx q[2];
rz(-2.0187536) q[2];
sx q[2];
rz(3.0222859) q[2];
rz(-1.9909158) q[3];
sx q[3];
rz(-2.2169866) q[3];
sx q[3];
rz(0.22056072) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
