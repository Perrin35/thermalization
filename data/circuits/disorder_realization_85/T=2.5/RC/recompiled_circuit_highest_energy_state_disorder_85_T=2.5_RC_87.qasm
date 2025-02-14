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
rz(0.97419089) q[0];
sx q[0];
rz(3.9904896) q[0];
sx q[0];
rz(11.150802) q[0];
rz(2.2805136) q[1];
sx q[1];
rz(-1.5134892) q[1];
sx q[1];
rz(2.3754062) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4199894) q[0];
sx q[0];
rz(-3.0675305) q[0];
sx q[0];
rz(-0.5038528) q[0];
x q[1];
rz(-3.1076961) q[2];
sx q[2];
rz(-1.4501415) q[2];
sx q[2];
rz(-1.1519037) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79957818) q[1];
sx q[1];
rz(-1.3675081) q[1];
sx q[1];
rz(-1.5631097) q[1];
x q[2];
rz(1.9835408) q[3];
sx q[3];
rz(-1.4719324) q[3];
sx q[3];
rz(-0.017798558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2506931) q[2];
sx q[2];
rz(-0.74960274) q[2];
sx q[2];
rz(3.0787943) q[2];
rz(1.7167669) q[3];
sx q[3];
rz(-1.3751605) q[3];
sx q[3];
rz(-0.73960251) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0733136) q[0];
sx q[0];
rz(-0.24389076) q[0];
sx q[0];
rz(-2.0252315) q[0];
rz(1.5694537) q[1];
sx q[1];
rz(-0.29466033) q[1];
sx q[1];
rz(1.1192628) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1946487) q[0];
sx q[0];
rz(-2.4341321) q[0];
sx q[0];
rz(-0.0012673541) q[0];
x q[1];
rz(-2.0541899) q[2];
sx q[2];
rz(-1.2366939) q[2];
sx q[2];
rz(-2.7887864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3685064) q[1];
sx q[1];
rz(-2.0952053) q[1];
sx q[1];
rz(1.6462516) q[1];
x q[2];
rz(2.4495115) q[3];
sx q[3];
rz(-1.0861703) q[3];
sx q[3];
rz(1.6177819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0548627) q[2];
sx q[2];
rz(-0.77069288) q[2];
sx q[2];
rz(1.9697624) q[2];
rz(2.150676) q[3];
sx q[3];
rz(-2.1358229) q[3];
sx q[3];
rz(1.2539554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1427926) q[0];
sx q[0];
rz(-0.87738377) q[0];
sx q[0];
rz(-3.0679833) q[0];
rz(2.7403846) q[1];
sx q[1];
rz(-2.7707272) q[1];
sx q[1];
rz(-2.4295095) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7382848) q[0];
sx q[0];
rz(-1.6845595) q[0];
sx q[0];
rz(3.0902181) q[0];
rz(-2.4855494) q[2];
sx q[2];
rz(-1.5719766) q[2];
sx q[2];
rz(-0.16108433) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5873196) q[1];
sx q[1];
rz(-1.463542) q[1];
sx q[1];
rz(2.8337906) q[1];
rz(1.6693657) q[3];
sx q[3];
rz(-0.53547066) q[3];
sx q[3];
rz(-0.55619538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.74505836) q[2];
sx q[2];
rz(-1.859954) q[2];
sx q[2];
rz(-2.7520531) q[2];
rz(2.7430096) q[3];
sx q[3];
rz(-1.7321209) q[3];
sx q[3];
rz(1.4220953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013537708) q[0];
sx q[0];
rz(-0.062352926) q[0];
sx q[0];
rz(2.7259977) q[0];
rz(1.2526814) q[1];
sx q[1];
rz(-2.0038192) q[1];
sx q[1];
rz(0.91928732) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18596953) q[0];
sx q[0];
rz(-2.9730995) q[0];
sx q[0];
rz(1.3244434) q[0];
rz(-pi) q[1];
rz(1.5933199) q[2];
sx q[2];
rz(-1.3078346) q[2];
sx q[2];
rz(2.7561273) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.63863149) q[1];
sx q[1];
rz(-2.161484) q[1];
sx q[1];
rz(0.69815127) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1197129) q[3];
sx q[3];
rz(-0.85441426) q[3];
sx q[3];
rz(2.7855887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.84819841) q[2];
sx q[2];
rz(-1.6323117) q[2];
sx q[2];
rz(-2.1634114) q[2];
rz(2.6350422) q[3];
sx q[3];
rz(-1.4975486) q[3];
sx q[3];
rz(-1.0165366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070234805) q[0];
sx q[0];
rz(-2.4429584) q[0];
sx q[0];
rz(0.51682669) q[0];
rz(1.0454073) q[1];
sx q[1];
rz(-2.1357336) q[1];
sx q[1];
rz(-0.4019092) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5346679) q[0];
sx q[0];
rz(-0.87550301) q[0];
sx q[0];
rz(0.57906998) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2847309) q[2];
sx q[2];
rz(-1.7334043) q[2];
sx q[2];
rz(2.3608077) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4337804) q[1];
sx q[1];
rz(-1.3715327) q[1];
sx q[1];
rz(0.0043745478) q[1];
rz(-pi) q[2];
rz(3.0203781) q[3];
sx q[3];
rz(-0.95887163) q[3];
sx q[3];
rz(2.2387089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4154633) q[2];
sx q[2];
rz(-1.9620506) q[2];
sx q[2];
rz(2.9187875) q[2];
rz(-0.86083096) q[3];
sx q[3];
rz(-1.6313044) q[3];
sx q[3];
rz(1.2739325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.05805) q[0];
sx q[0];
rz(-1.2419751) q[0];
sx q[0];
rz(0.85839957) q[0];
rz(-2.3421613) q[1];
sx q[1];
rz(-2.1245427) q[1];
sx q[1];
rz(1.7564836) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1027062) q[0];
sx q[0];
rz(-1.2974076) q[0];
sx q[0];
rz(-2.9486743) q[0];
rz(-pi) q[1];
rz(-3.0869003) q[2];
sx q[2];
rz(-1.8963976) q[2];
sx q[2];
rz(-0.97132746) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0771602) q[1];
sx q[1];
rz(-0.61187498) q[1];
sx q[1];
rz(2.9303657) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0812182) q[3];
sx q[3];
rz(-1.7822232) q[3];
sx q[3];
rz(-1.2871773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.73056209) q[2];
sx q[2];
rz(-2.4081814) q[2];
sx q[2];
rz(2.9900271) q[2];
rz(-1.1847121) q[3];
sx q[3];
rz(-2.0047174) q[3];
sx q[3];
rz(-2.241551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2043692) q[0];
sx q[0];
rz(-1.3441514) q[0];
sx q[0];
rz(-0.7835266) q[0];
rz(-2.9339058) q[1];
sx q[1];
rz(-0.91778225) q[1];
sx q[1];
rz(-1.2088561) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1087813) q[0];
sx q[0];
rz(-2.0447123) q[0];
sx q[0];
rz(1.117871) q[0];
rz(-pi) q[1];
rz(0.54401347) q[2];
sx q[2];
rz(-0.90780963) q[2];
sx q[2];
rz(1.3617409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5294339) q[1];
sx q[1];
rz(-0.72996695) q[1];
sx q[1];
rz(1.5117701) q[1];
rz(-1.1829258) q[3];
sx q[3];
rz(-0.3976477) q[3];
sx q[3];
rz(-2.7499547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2850767) q[2];
sx q[2];
rz(-2.3975942) q[2];
sx q[2];
rz(3.0169955) q[2];
rz(1.3397269) q[3];
sx q[3];
rz(-2.7281269) q[3];
sx q[3];
rz(-1.193803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0143921) q[0];
sx q[0];
rz(-1.7008282) q[0];
sx q[0];
rz(0.097323962) q[0];
rz(-0.65517455) q[1];
sx q[1];
rz(-0.56071463) q[1];
sx q[1];
rz(-2.4526144) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45132056) q[0];
sx q[0];
rz(-0.49271944) q[0];
sx q[0];
rz(-1.7283597) q[0];
rz(-pi) q[1];
rz(-1.6332878) q[2];
sx q[2];
rz(-0.35258503) q[2];
sx q[2];
rz(3.0931713) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0216768) q[1];
sx q[1];
rz(-2.4525453) q[1];
sx q[1];
rz(0.20153789) q[1];
rz(0.68920691) q[3];
sx q[3];
rz(-1.4446994) q[3];
sx q[3];
rz(1.4339303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24136209) q[2];
sx q[2];
rz(-2.1817544) q[2];
sx q[2];
rz(2.0755365) q[2];
rz(0.42029941) q[3];
sx q[3];
rz(-1.2371141) q[3];
sx q[3];
rz(1.458781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083753839) q[0];
sx q[0];
rz(-0.76863113) q[0];
sx q[0];
rz(0.88430697) q[0];
rz(-0.090944313) q[1];
sx q[1];
rz(-2.2021286) q[1];
sx q[1];
rz(-1.1788064) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1183848) q[0];
sx q[0];
rz(-1.7309446) q[0];
sx q[0];
rz(1.5051756) q[0];
x q[1];
rz(-0.0021807533) q[2];
sx q[2];
rz(-2.2366887) q[2];
sx q[2];
rz(-1.8357753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2388666) q[1];
sx q[1];
rz(-1.7416493) q[1];
sx q[1];
rz(2.577223) q[1];
rz(-1.6903773) q[3];
sx q[3];
rz(-1.371583) q[3];
sx q[3];
rz(-1.3026893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74943298) q[2];
sx q[2];
rz(-1.1460816) q[2];
sx q[2];
rz(2.667099) q[2];
rz(0.042304603) q[3];
sx q[3];
rz(-1.5174805) q[3];
sx q[3];
rz(1.8208985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.6831191) q[0];
sx q[0];
rz(-1.1011769) q[0];
sx q[0];
rz(1.3018357) q[0];
rz(2.5104972) q[1];
sx q[1];
rz(-1.0429691) q[1];
sx q[1];
rz(1.016681) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23970579) q[0];
sx q[0];
rz(-1.8586057) q[0];
sx q[0];
rz(-2.2255656) q[0];
rz(-pi) q[1];
x q[1];
rz(0.053875845) q[2];
sx q[2];
rz(-1.6695256) q[2];
sx q[2];
rz(1.7629262) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70003033) q[1];
sx q[1];
rz(-2.4228281) q[1];
sx q[1];
rz(1.5573386) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74211699) q[3];
sx q[3];
rz(-0.85608965) q[3];
sx q[3];
rz(-1.2679188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8299228) q[2];
sx q[2];
rz(-2.416478) q[2];
sx q[2];
rz(-0.97490087) q[2];
rz(1.7389065) q[3];
sx q[3];
rz(-2.1699173) q[3];
sx q[3];
rz(-1.5909125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.0065895157) q[0];
sx q[0];
rz(-2.4000744) q[0];
sx q[0];
rz(0.97493521) q[0];
rz(1.4047752) q[1];
sx q[1];
rz(-1.4092696) q[1];
sx q[1];
rz(-0.41913941) q[1];
rz(-2.1640833) q[2];
sx q[2];
rz(-2.0089663) q[2];
sx q[2];
rz(0.51539863) q[2];
rz(2.9007965) q[3];
sx q[3];
rz(-0.44178648) q[3];
sx q[3];
rz(1.0159258) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
