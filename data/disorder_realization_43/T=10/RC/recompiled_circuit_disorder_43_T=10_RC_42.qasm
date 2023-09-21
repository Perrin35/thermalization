OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(1.8068846) q[0];
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(0.97595739) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1057518) q[0];
sx q[0];
rz(-2.5076206) q[0];
sx q[0];
rz(-0.41497725) q[0];
x q[1];
rz(-1.0693597) q[2];
sx q[2];
rz(-0.62383365) q[2];
sx q[2];
rz(1.8471579) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18260278) q[1];
sx q[1];
rz(-2.5187153) q[1];
sx q[1];
rz(1.3401003) q[1];
x q[2];
rz(-2.3372041) q[3];
sx q[3];
rz(-1.2803004) q[3];
sx q[3];
rz(-0.94798541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7261937) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(2.3577918) q[2];
rz(-0.33256724) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6587104) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(2.9630307) q[0];
rz(-1.3372955) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(3.1352502) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7550678) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(0.13271876) q[0];
rz(2.6132934) q[2];
sx q[2];
rz(-1.0353147) q[2];
sx q[2];
rz(-2.2250125) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.43286846) q[1];
sx q[1];
rz(-1.2607062) q[1];
sx q[1];
rz(1.5061782) q[1];
rz(-pi) q[2];
rz(-2.5793521) q[3];
sx q[3];
rz(-1.3573109) q[3];
sx q[3];
rz(-1.7052887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94770849) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-0.87810278) q[2];
rz(-2.2795423) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(-0.0058962065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55643117) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(0.01097824) q[0];
rz(2.7745461) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(3.045851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84155267) q[0];
sx q[0];
rz(-1.3944355) q[0];
sx q[0];
rz(-0.15533133) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0289621) q[2];
sx q[2];
rz(-2.0320227) q[2];
sx q[2];
rz(-1.065965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9853471) q[1];
sx q[1];
rz(-0.74838446) q[1];
sx q[1];
rz(-2.2845539) q[1];
rz(-pi) q[2];
rz(2.5903969) q[3];
sx q[3];
rz(-1.5407729) q[3];
sx q[3];
rz(2.4195645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0720955) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(0.83479184) q[2];
rz(0.21162027) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(0.37477469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(-2.7950177) q[0];
rz(-2.6158781) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(-2.1077164) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047065145) q[0];
sx q[0];
rz(-2.3405502) q[0];
sx q[0];
rz(0.69181504) q[0];
rz(0.526555) q[2];
sx q[2];
rz(-1.0609846) q[2];
sx q[2];
rz(-1.1916135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.911072) q[1];
sx q[1];
rz(-1.6358747) q[1];
sx q[1];
rz(-1.3793524) q[1];
rz(2.8598966) q[3];
sx q[3];
rz(-1.8857737) q[3];
sx q[3];
rz(1.7236934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6870849) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.7948077) q[2];
rz(0.421031) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(-0.69798654) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43276697) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(-0.41473266) q[0];
rz(-1.3955836) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(2.5674852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2643471) q[0];
sx q[0];
rz(-2.1640722) q[0];
sx q[0];
rz(-0.1546774) q[0];
rz(1.5933883) q[2];
sx q[2];
rz(-2.2859757) q[2];
sx q[2];
rz(1.9695645) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2327323) q[1];
sx q[1];
rz(-2.6870011) q[1];
sx q[1];
rz(-2.4652387) q[1];
x q[2];
rz(2.2361034) q[3];
sx q[3];
rz(-1.724616) q[3];
sx q[3];
rz(-0.32115768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2888912) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(1.627702) q[2];
rz(3.1001575) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(-0.07853011) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067327499) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(2.4434027) q[0];
rz(-0.16695887) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(2.4093157) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9376611) q[0];
sx q[0];
rz(-1.7294149) q[0];
sx q[0];
rz(-0.21477867) q[0];
rz(-pi) q[1];
rz(-0.68533021) q[2];
sx q[2];
rz(-0.88980674) q[2];
sx q[2];
rz(0.53390098) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4477168) q[1];
sx q[1];
rz(-0.69732053) q[1];
sx q[1];
rz(-2.7867774) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5712588) q[3];
sx q[3];
rz(-1.3864281) q[3];
sx q[3];
rz(2.1031228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7531062) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(1.1614655) q[2];
rz(-0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(-1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5724065) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(-2.9851595) q[0];
rz(-2.6898443) q[1];
sx q[1];
rz(-0.86507559) q[1];
sx q[1];
rz(-2.3715473) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55243385) q[0];
sx q[0];
rz(-1.7434412) q[0];
sx q[0];
rz(0.60683672) q[0];
rz(-pi) q[1];
rz(0.40600834) q[2];
sx q[2];
rz(-1.0205262) q[2];
sx q[2];
rz(1.5232435) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5321977) q[1];
sx q[1];
rz(-0.97071338) q[1];
sx q[1];
rz(-2.1983912) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2392063) q[3];
sx q[3];
rz(-1.0911897) q[3];
sx q[3];
rz(2.1945206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4975171) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(1.0151803) q[2];
rz(1.7049568) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(-0.59593433) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5148233) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(0.24169895) q[0];
rz(-0.73879755) q[1];
sx q[1];
rz(-2.6530478) q[1];
sx q[1];
rz(-0.36639211) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9241656) q[0];
sx q[0];
rz(-2.1806742) q[0];
sx q[0];
rz(0.92745552) q[0];
rz(-pi) q[1];
rz(1.1364469) q[2];
sx q[2];
rz(-2.1546116) q[2];
sx q[2];
rz(-2.7455612) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.061325039) q[1];
sx q[1];
rz(-0.72239164) q[1];
sx q[1];
rz(2.8894043) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4922769) q[3];
sx q[3];
rz(-0.55320569) q[3];
sx q[3];
rz(-1.6012524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.9423449) q[2];
sx q[2];
rz(-2.8477342) q[2];
rz(0.014523225) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83207399) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(-2.2576387) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(0.79137897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6521097) q[0];
sx q[0];
rz(-1.2533422) q[0];
sx q[0];
rz(-0.86274685) q[0];
rz(1.0857401) q[2];
sx q[2];
rz(-1.1527449) q[2];
sx q[2];
rz(0.56318356) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1014175) q[1];
sx q[1];
rz(-1.1244332) q[1];
sx q[1];
rz(1.1091713) q[1];
x q[2];
rz(1.8459686) q[3];
sx q[3];
rz(-2.5857946) q[3];
sx q[3];
rz(-2.9142227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.066594921) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(-2.1203314) q[2];
rz(-0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(-0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1148949) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(0.60780203) q[0];
rz(-0.23884493) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(-1.3806608) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62771195) q[0];
sx q[0];
rz(-2.4988334) q[0];
sx q[0];
rz(1.7895262) q[0];
x q[1];
rz(-0.025823921) q[2];
sx q[2];
rz(-1.1354453) q[2];
sx q[2];
rz(-0.045407427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3874515) q[1];
sx q[1];
rz(-1.2536465) q[1];
sx q[1];
rz(-1.1412568) q[1];
rz(0.72762604) q[3];
sx q[3];
rz(-0.25087038) q[3];
sx q[3];
rz(-2.6550456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7878788) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-0.33622462) q[2];
rz(2.0119038) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1375785) q[0];
sx q[0];
rz(-1.4083569) q[0];
sx q[0];
rz(1.2271723) q[0];
rz(2.5972988) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(-0.82472807) q[2];
sx q[2];
rz(-1.8246973) q[2];
sx q[2];
rz(1.4551103) q[2];
rz(-1.0014793) q[3];
sx q[3];
rz(-1.7703198) q[3];
sx q[3];
rz(1.9386335) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];