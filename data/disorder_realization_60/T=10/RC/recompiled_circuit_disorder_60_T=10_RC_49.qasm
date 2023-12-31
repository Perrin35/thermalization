OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1922167) q[0];
sx q[0];
rz(-2.0944216) q[0];
sx q[0];
rz(3.0728683) q[0];
rz(1.7460495) q[1];
sx q[1];
rz(-1.6092602) q[1];
sx q[1];
rz(-1.2083763) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87586227) q[0];
sx q[0];
rz(-3.0348572) q[0];
sx q[0];
rz(2.021832) q[0];
rz(-pi) q[1];
rz(3/(10*pi)) q[2];
sx q[2];
rz(-0.12636939) q[2];
sx q[2];
rz(0.40590826) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0304058) q[1];
sx q[1];
rz(-2.5136247) q[1];
sx q[1];
rz(-0.18917947) q[1];
rz(1.0923907) q[3];
sx q[3];
rz(-2.0432825) q[3];
sx q[3];
rz(2.9977968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8157114) q[2];
sx q[2];
rz(-1.7408966) q[2];
sx q[2];
rz(1.4665843) q[2];
rz(0.69774929) q[3];
sx q[3];
rz(-1.1013228) q[3];
sx q[3];
rz(2.3944323) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0242457) q[0];
sx q[0];
rz(-1.1276561) q[0];
sx q[0];
rz(1.9673989) q[0];
rz(-0.17114561) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(2.8443964) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8288119) q[0];
sx q[0];
rz(-0.096185616) q[0];
sx q[0];
rz(1.0440774) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0969639) q[2];
sx q[2];
rz(-2.2261438) q[2];
sx q[2];
rz(1.6608134) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.17644037) q[1];
sx q[1];
rz(-2.0808176) q[1];
sx q[1];
rz(-2.7494207) q[1];
rz(0.76241775) q[3];
sx q[3];
rz(-1.276351) q[3];
sx q[3];
rz(0.03014119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.16195665) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(0.61398181) q[2];
rz(0.87614122) q[3];
sx q[3];
rz(-2.4738779) q[3];
sx q[3];
rz(-2.5045625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0959594) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(-2.8339548) q[0];
rz(-2.3930507) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(0.83980733) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52536406) q[0];
sx q[0];
rz(-1.9814241) q[0];
sx q[0];
rz(-1.8812268) q[0];
rz(-pi) q[1];
rz(0.26250458) q[2];
sx q[2];
rz(-2.7911107) q[2];
sx q[2];
rz(-1.0002491) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0963904) q[1];
sx q[1];
rz(-0.55711105) q[1];
sx q[1];
rz(0.60369173) q[1];
rz(-pi) q[2];
rz(2.2037376) q[3];
sx q[3];
rz(-1.4301392) q[3];
sx q[3];
rz(-1.3682287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.86747375) q[2];
sx q[2];
rz(-1.7904736) q[2];
sx q[2];
rz(-1.8236558) q[2];
rz(1.9258202) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(1.5095476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3777305) q[0];
sx q[0];
rz(-1.7611935) q[0];
sx q[0];
rz(-0.44556251) q[0];
rz(2.5207649) q[1];
sx q[1];
rz(-1.2460243) q[1];
sx q[1];
rz(-2.1760118) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5774028) q[0];
sx q[0];
rz(-1.8468542) q[0];
sx q[0];
rz(0.77703707) q[0];
rz(-1.6057274) q[2];
sx q[2];
rz(-0.72872773) q[2];
sx q[2];
rz(-2.5503089) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1196739) q[1];
sx q[1];
rz(-1.5389171) q[1];
sx q[1];
rz(1.8244513) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86430092) q[3];
sx q[3];
rz(-0.14926499) q[3];
sx q[3];
rz(2.5316558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.821637) q[2];
sx q[2];
rz(-2.379202) q[2];
sx q[2];
rz(-2.2122673) q[2];
rz(2.4980513) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(-2.1381366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-1.6699566) q[0];
sx q[0];
rz(-1.5595373) q[0];
sx q[0];
rz(1.2840282) q[0];
rz(0.28981003) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(1.0481542) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5149887) q[0];
sx q[0];
rz(-2.0467313) q[0];
sx q[0];
rz(2.3355052) q[0];
rz(-pi) q[1];
rz(-0.83277793) q[2];
sx q[2];
rz(-1.1156429) q[2];
sx q[2];
rz(0.27311329) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3874665) q[1];
sx q[1];
rz(-1.4782463) q[1];
sx q[1];
rz(-2.8184163) q[1];
rz(0.48091472) q[3];
sx q[3];
rz(-0.9489343) q[3];
sx q[3];
rz(2.8180168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2709048) q[2];
sx q[2];
rz(-1.015816) q[2];
sx q[2];
rz(-2.2407545) q[2];
rz(1.0926931) q[3];
sx q[3];
rz(-1.0034424) q[3];
sx q[3];
rz(1.9074915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822405) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(-0.59610468) q[0];
rz(-1.4959363) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(-1.8966819) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0672452) q[0];
sx q[0];
rz(-1.0676358) q[0];
sx q[0];
rz(-2.7748681) q[0];
rz(1.2539358) q[2];
sx q[2];
rz(-1.1057901) q[2];
sx q[2];
rz(2.3777865) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5925496) q[1];
sx q[1];
rz(-0.92714308) q[1];
sx q[1];
rz(-2.2085269) q[1];
rz(-pi) q[2];
rz(-0.2480898) q[3];
sx q[3];
rz(-1.7956942) q[3];
sx q[3];
rz(3.0450862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.231455) q[2];
sx q[2];
rz(-0.69880501) q[2];
sx q[2];
rz(-0.22496741) q[2];
rz(-3.0531626) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(-2.6627873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.6769619) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(0.27994573) q[0];
rz(1.4631368) q[1];
sx q[1];
rz(-1.2660374) q[1];
sx q[1];
rz(2.8889012) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42851105) q[0];
sx q[0];
rz(-1.6185074) q[0];
sx q[0];
rz(-1.2776572) q[0];
x q[1];
rz(-1.5451317) q[2];
sx q[2];
rz(-0.63458323) q[2];
sx q[2];
rz(-2.6013825) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.686324) q[1];
sx q[1];
rz(-0.9372006) q[1];
sx q[1];
rz(0.92647657) q[1];
rz(-pi) q[2];
rz(-2.2884376) q[3];
sx q[3];
rz(-0.95082885) q[3];
sx q[3];
rz(2.6593047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32020405) q[2];
sx q[2];
rz(-0.89367047) q[2];
sx q[2];
rz(-1.6097216) q[2];
rz(1.1931217) q[3];
sx q[3];
rz(-0.90819287) q[3];
sx q[3];
rz(-1.0866722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10483345) q[0];
sx q[0];
rz(-1.7585254) q[0];
sx q[0];
rz(1.4861134) q[0];
rz(0.2688109) q[1];
sx q[1];
rz(-2.0116282) q[1];
sx q[1];
rz(-0.2789467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7689965) q[0];
sx q[0];
rz(-1.0770174) q[0];
sx q[0];
rz(-0.84817024) q[0];
rz(-1.0686764) q[2];
sx q[2];
rz(-1.7218105) q[2];
sx q[2];
rz(2.6627024) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8262417) q[1];
sx q[1];
rz(-0.8289753) q[1];
sx q[1];
rz(0.51893236) q[1];
x q[2];
rz(1.7989743) q[3];
sx q[3];
rz(-0.5842714) q[3];
sx q[3];
rz(-2.1669471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3387317) q[2];
sx q[2];
rz(-1.4638476) q[2];
sx q[2];
rz(-1.5926682) q[2];
rz(-1.0507978) q[3];
sx q[3];
rz(-1.2069586) q[3];
sx q[3];
rz(1.1999493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779697) q[0];
sx q[0];
rz(-2.2080053) q[0];
sx q[0];
rz(1.2532225) q[0];
rz(-0.62250096) q[1];
sx q[1];
rz(-1.4600735) q[1];
sx q[1];
rz(-1.1463096) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2749467) q[0];
sx q[0];
rz(-1.5005611) q[0];
sx q[0];
rz(-2.3932891) q[0];
x q[1];
rz(-2.3290645) q[2];
sx q[2];
rz(-2.3000237) q[2];
sx q[2];
rz(1.3844045) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6282564) q[1];
sx q[1];
rz(-2.4504821) q[1];
sx q[1];
rz(-1.4555132) q[1];
rz(-1.1504437) q[3];
sx q[3];
rz(-0.97878362) q[3];
sx q[3];
rz(-1.0279442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.15381947) q[2];
sx q[2];
rz(-1.035707) q[2];
sx q[2];
rz(2.1257607) q[2];
rz(2.231797) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(-0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0125473) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(-3.1273499) q[0];
rz(-0.8447389) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(-1.7600118) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81603564) q[0];
sx q[0];
rz(-2.9986458) q[0];
sx q[0];
rz(-1.4965034) q[0];
x q[1];
rz(-3.0773452) q[2];
sx q[2];
rz(-0.99284961) q[2];
sx q[2];
rz(-0.77677514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6784918) q[1];
sx q[1];
rz(-1.4271724) q[1];
sx q[1];
rz(-2.1302057) q[1];
rz(-pi) q[2];
rz(0.53346177) q[3];
sx q[3];
rz(-1.3812314) q[3];
sx q[3];
rz(1.61434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0766729) q[2];
sx q[2];
rz(-1.939247) q[2];
sx q[2];
rz(0.98999611) q[2];
rz(2.2475217) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(0.22542424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0902696) q[0];
sx q[0];
rz(-2.4933503) q[0];
sx q[0];
rz(-1.1052263) q[0];
rz(-1.3399711) q[1];
sx q[1];
rz(-2.5201288) q[1];
sx q[1];
rz(-2.7531243) q[1];
rz(0.71628911) q[2];
sx q[2];
rz(-1.204797) q[2];
sx q[2];
rz(2.7439678) q[2];
rz(2.1918478) q[3];
sx q[3];
rz(-1.2093778) q[3];
sx q[3];
rz(-2.9686684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
