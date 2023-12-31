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
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(-2.1656353) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0358409) q[0];
sx q[0];
rz(-2.5076206) q[0];
sx q[0];
rz(0.41497725) q[0];
rz(-0.33306723) q[2];
sx q[2];
rz(-2.1085848) q[2];
sx q[2];
rz(1.8884459) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9589899) q[1];
sx q[1];
rz(-2.5187153) q[1];
sx q[1];
rz(1.3401003) q[1];
rz(-0.39335143) q[3];
sx q[3];
rz(-2.297612) q[3];
sx q[3];
rz(-2.2497183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(0.78380084) q[2];
rz(2.8090254) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6587104) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(2.9630307) q[0];
rz(-1.8042971) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(3.1352502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7550678) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(-3.0088739) q[0];
rz(-pi) q[1];
rz(0.86654051) q[2];
sx q[2];
rz(-0.73359493) q[2];
sx q[2];
rz(-1.7689592) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22390631) q[1];
sx q[1];
rz(-0.31653857) q[1];
sx q[1];
rz(2.9427337) q[1];
rz(2.7553495) q[3];
sx q[3];
rz(-0.59730232) q[3];
sx q[3];
rz(2.951705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1938842) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(2.2795423) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(-0.0058962065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.55643117) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(0.01097824) q[0];
rz(2.7745461) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(-0.095741622) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11287963) q[0];
sx q[0];
rz(-2.9071147) q[0];
sx q[0];
rz(-2.2857091) q[0];
x q[1];
rz(-0.7272561) q[2];
sx q[2];
rz(-0.63823344) q[2];
sx q[2];
rz(-2.9122695) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.84903753) q[1];
sx q[1];
rz(-2.0325066) q[1];
sx q[1];
rz(0.95878102) q[1];
x q[2];
rz(-3.0843094) q[3];
sx q[3];
rz(-2.5896642) q[3];
sx q[3];
rz(-2.3416167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0694971) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(0.83479184) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(0.37477469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.19462207) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(2.7950177) q[0];
rz(2.6158781) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(2.1077164) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3165986) q[0];
sx q[0];
rz(-0.98485095) q[0];
sx q[0];
rz(-0.98866776) q[0];
rz(-pi) q[1];
rz(2.6150377) q[2];
sx q[2];
rz(-1.0609846) q[2];
sx q[2];
rz(1.1916135) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8139207) q[1];
sx q[1];
rz(-1.76183) q[1];
sx q[1];
rz(-0.066285985) q[1];
rz(-0.86446188) q[3];
sx q[3];
rz(-0.41941386) q[3];
sx q[3];
rz(2.1692587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6870849) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.7948077) q[2];
rz(-0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(2.4436061) q[3];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088257) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(0.41473266) q[0];
rz(-1.3955836) q[1];
sx q[1];
rz(-0.65892977) q[1];
sx q[1];
rz(-2.5674852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78050437) q[0];
sx q[0];
rz(-1.6988806) q[0];
sx q[0];
rz(2.1696521) q[0];
x q[1];
rz(0.026002361) q[2];
sx q[2];
rz(-2.4261195) q[2];
sx q[2];
rz(-2.0040087) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.50358665) q[1];
sx q[1];
rz(-1.2212911) q[1];
sx q[1];
rz(1.2739146) q[1];
rz(-0.19458171) q[3];
sx q[3];
rz(-0.91472018) q[3];
sx q[3];
rz(1.7723099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2888912) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(-1.627702) q[2];
rz(0.04143516) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(0.07853011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067327499) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(-0.69818991) q[0];
rz(-0.16695887) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(-0.73227698) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9376611) q[0];
sx q[0];
rz(-1.4121778) q[0];
sx q[0];
rz(0.21477867) q[0];
x q[1];
rz(0.68533021) q[2];
sx q[2];
rz(-0.88980674) q[2];
sx q[2];
rz(-0.53390098) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4477168) q[1];
sx q[1];
rz(-2.4442721) q[1];
sx q[1];
rz(-0.35481528) q[1];
rz(1.7888277) q[3];
sx q[3];
rz(-2.1302967) q[3];
sx q[3];
rz(-2.7263209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3884864) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(1.1614655) q[2];
rz(2.8325864) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.5724065) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(-0.15643315) q[0];
rz(-2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(2.3715473) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8806481) q[0];
sx q[0];
rz(-0.62793193) q[0];
sx q[0];
rz(2.8448366) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7355843) q[2];
sx q[2];
rz(-1.0205262) q[2];
sx q[2];
rz(-1.5232435) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7913831) q[1];
sx q[1];
rz(-2.0767127) q[1];
sx q[1];
rz(-2.4398068) q[1];
rz(-1.2392063) q[3];
sx q[3];
rz(-1.0911897) q[3];
sx q[3];
rz(-0.94707205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4975171) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(1.0151803) q[2];
rz(-1.4366359) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.5148233) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(0.24169895) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-0.36639211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9241656) q[0];
sx q[0];
rz(-0.96091849) q[0];
sx q[0];
rz(-0.92745552) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5744152) q[2];
sx q[2];
rz(-2.4293373) q[2];
sx q[2];
rz(2.8380307) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3185127) q[1];
sx q[1];
rz(-1.4050583) q[1];
sx q[1];
rz(-2.4351099) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6493158) q[3];
sx q[3];
rz(-2.588387) q[3];
sx q[3];
rz(-1.6012524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1239132) q[2];
sx q[2];
rz(-1.9423449) q[2];
sx q[2];
rz(-2.8477342) q[2];
rz(-3.1270694) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83207399) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(0.88395399) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(-0.79137897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26830772) q[0];
sx q[0];
rz(-0.76457667) q[0];
sx q[0];
rz(1.1029878) q[0];
rz(0.80957885) q[2];
sx q[2];
rz(-0.62925816) q[2];
sx q[2];
rz(1.4776243) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1014175) q[1];
sx q[1];
rz(-1.1244332) q[1];
sx q[1];
rz(2.0324213) q[1];
rz(-pi) q[2];
rz(1.295624) q[3];
sx q[3];
rz(-0.55579805) q[3];
sx q[3];
rz(-2.9142227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0749977) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(-2.1203314) q[2];
rz(-0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(2.5436201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1148949) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(0.60780203) q[0];
rz(-0.23884493) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(-1.3806608) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5138807) q[0];
sx q[0];
rz(-2.4988334) q[0];
sx q[0];
rz(1.7895262) q[0];
rz(-pi) q[1];
rz(1.1353178) q[2];
sx q[2];
rz(-1.594211) q[2];
sx q[2];
rz(-1.6053111) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32523649) q[1];
sx q[1];
rz(-1.9776077) q[1];
sx q[1];
rz(-0.34646323) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72762604) q[3];
sx q[3];
rz(-2.8907223) q[3];
sx q[3];
rz(2.6550456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7878788) q[2];
sx q[2];
rz(-1.6106662) q[2];
sx q[2];
rz(-2.805368) q[2];
rz(1.1296889) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(-2.0000134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1375785) q[0];
sx q[0];
rz(-1.4083569) q[0];
sx q[0];
rz(1.2271723) q[0];
rz(0.54429383) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(-2.8019194) q[2];
sx q[2];
rz(-0.85396955) q[2];
sx q[2];
rz(0.11243482) q[2];
rz(-1.2119157) q[3];
sx q[3];
rz(-2.5419895) q[3];
sx q[3];
rz(0.06751577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
