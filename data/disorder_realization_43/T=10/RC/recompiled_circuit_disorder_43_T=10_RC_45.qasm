OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2551978) q[0];
sx q[0];
rz(4.3918443) q[0];
sx q[0];
rz(10.759486) q[0];
rz(2.788738) q[1];
sx q[1];
rz(3.3021441) q[1];
sx q[1];
rz(8.4488206) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87603509) q[0];
sx q[0];
rz(-1.3296488) q[0];
sx q[0];
rz(-2.5493456) q[0];
rz(-pi) q[1];
rz(1.0693597) q[2];
sx q[2];
rz(-0.62383365) q[2];
sx q[2];
rz(-1.8471579) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.576697) q[1];
sx q[1];
rz(-1.4370059) q[1];
sx q[1];
rz(0.96058515) q[1];
rz(1.9777771) q[3];
sx q[3];
rz(-2.3325936) q[3];
sx q[3];
rz(0.33363261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41539899) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(-0.78380084) q[2];
rz(2.8090254) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587104) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(-0.17856199) q[0];
rz(1.3372955) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(-3.1352502) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7550678) q[0];
sx q[0];
rz(-1.3773943) q[0];
sx q[0];
rz(0.13271876) q[0];
x q[1];
rz(2.2750521) q[2];
sx q[2];
rz(-2.4079977) q[2];
sx q[2];
rz(1.3726335) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22390631) q[1];
sx q[1];
rz(-2.8250541) q[1];
sx q[1];
rz(0.19885893) q[1];
rz(1.3199602) q[3];
sx q[3];
rz(-1.0227961) q[3];
sx q[3];
rz(2.8743924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.94770849) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(-2.2795423) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(0.0058962065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55643117) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(0.01097824) q[0];
rz(-0.36704656) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(3.045851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.028713) q[0];
sx q[0];
rz(-0.23447795) q[0];
sx q[0];
rz(2.2857091) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7272561) q[2];
sx q[2];
rz(-2.5033592) q[2];
sx q[2];
rz(-0.22932316) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15624554) q[1];
sx q[1];
rz(-2.3932082) q[1];
sx q[1];
rz(2.2845539) q[1];
x q[2];
rz(0.5511958) q[3];
sx q[3];
rz(-1.6008198) q[3];
sx q[3];
rz(-0.72202819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0694971) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(2.3068008) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(-0.37477469) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(-0.34657493) q[0];
rz(-2.6158781) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(2.1077164) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047065145) q[0];
sx q[0];
rz(-0.80104242) q[0];
sx q[0];
rz(-0.69181504) q[0];
x q[1];
rz(-2.1448574) q[2];
sx q[2];
rz(-1.1166995) q[2];
sx q[2];
rz(0.10276375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8139207) q[1];
sx q[1];
rz(-1.76183) q[1];
sx q[1];
rz(3.0753067) q[1];
x q[2];
rz(1.8978118) q[3];
sx q[3];
rz(-1.3033086) q[3];
sx q[3];
rz(-0.063484065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.45450777) q[2];
sx q[2];
rz(-1.4322174) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(-0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(-0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43276697) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(-2.72686) q[0];
rz(1.3955836) q[1];
sx q[1];
rz(-0.65892977) q[1];
sx q[1];
rz(2.5674852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5363279) q[0];
sx q[0];
rz(-0.61075532) q[0];
sx q[0];
rz(1.7954134) q[0];
rz(-pi) q[1];
x q[1];
rz(0.026002361) q[2];
sx q[2];
rz(-2.4261195) q[2];
sx q[2];
rz(-2.0040087) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9628323) q[1];
sx q[1];
rz(-1.292359) q[1];
sx q[1];
rz(2.7774485) q[1];
rz(1.3247213) q[3];
sx q[3];
rz(-2.4613791) q[3];
sx q[3];
rz(-1.0567997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2888912) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(1.5138907) q[2];
rz(-3.1001575) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(0.07853011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067327499) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(0.69818991) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(-0.73227698) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0064272881) q[0];
sx q[0];
rz(-2.8753202) q[0];
sx q[0];
rz(-2.4977495) q[0];
rz(0.76262577) q[2];
sx q[2];
rz(-1.0566933) q[2];
sx q[2];
rz(-2.5800173) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.69387586) q[1];
sx q[1];
rz(-0.69732053) q[1];
sx q[1];
rz(-0.35481528) q[1];
rz(-2.5712588) q[3];
sx q[3];
rz(-1.7551646) q[3];
sx q[3];
rz(1.0384699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3884864) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(1.9801271) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56918615) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(-2.9851595) q[0];
rz(2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(-2.3715473) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89966398) q[0];
sx q[0];
rz(-0.97424346) q[0];
sx q[0];
rz(-1.3616189) q[0];
x q[1];
rz(0.98203512) q[2];
sx q[2];
rz(-1.9141478) q[2];
sx q[2];
rz(-2.8729168) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.60939497) q[1];
sx q[1];
rz(-2.1708793) q[1];
sx q[1];
rz(-0.9432015) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6386912) q[3];
sx q[3];
rz(-1.8637878) q[3];
sx q[3];
rz(-2.6754232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4975171) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(2.1264123) q[2];
rz(-1.4366359) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(-0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(2.7752005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0893433) q[0];
sx q[0];
rz(-1.0567259) q[0];
sx q[0];
rz(2.423717) q[0];
x q[1];
rz(2.5744152) q[2];
sx q[2];
rz(-2.4293373) q[2];
sx q[2];
rz(-2.8380307) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7494292) q[1];
sx q[1];
rz(-0.87595075) q[1];
sx q[1];
rz(-1.3543345) q[1];
rz(-2.1226235) q[3];
sx q[3];
rz(-1.5295715) q[3];
sx q[3];
rz(0.036389694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(2.8477342) q[2];
rz(0.014523225) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(-0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095187) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(-2.2576387) q[0];
rz(2.4813095) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(-2.3502137) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.489483) q[0];
sx q[0];
rz(-1.2533422) q[0];
sx q[0];
rz(-0.86274685) q[0];
rz(1.0857401) q[2];
sx q[2];
rz(-1.1527449) q[2];
sx q[2];
rz(-2.5784091) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3190805) q[1];
sx q[1];
rz(-1.1573536) q[1];
sx q[1];
rz(0.49088571) q[1];
rz(-pi) q[2];
rz(-2.1095554) q[3];
sx q[3];
rz(-1.7146535) q[3];
sx q[3];
rz(-2.0335576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.066594921) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(1.0212612) q[2];
rz(0.3785454) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(-0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026697712) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(-0.60780203) q[0];
rz(0.23884493) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(1.7609319) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5138807) q[0];
sx q[0];
rz(-2.4988334) q[0];
sx q[0];
rz(1.3520665) q[0];
rz(-pi) q[1];
rz(3.1157687) q[2];
sx q[2];
rz(-2.0061473) q[2];
sx q[2];
rz(0.045407427) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32523649) q[1];
sx q[1];
rz(-1.9776077) q[1];
sx q[1];
rz(0.34646323) q[1];
rz(-pi) q[2];
rz(1.7396183) q[3];
sx q[3];
rz(-1.7572548) q[3];
sx q[3];
rz(2.8846915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3537139) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-2.805368) q[2];
rz(-1.1296889) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(2.0000134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0040141) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(2.5972988) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(1.9359679) q[2];
sx q[2];
rz(-2.361459) q[2];
sx q[2];
rz(-0.38103719) q[2];
rz(-0.23562283) q[3];
sx q[3];
rz(-2.1274673) q[3];
sx q[3];
rz(-2.6475788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
