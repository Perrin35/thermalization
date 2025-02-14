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
rz(0.089601547) q[0];
sx q[0];
rz(-0.90553415) q[0];
sx q[0];
rz(-0.18984689) q[0];
rz(2.7449961) q[1];
sx q[1];
rz(-0.34062579) q[1];
sx q[1];
rz(-0.61000282) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1413795) q[0];
sx q[0];
rz(-1.3249517) q[0];
sx q[0];
rz(2.8971998) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9556263) q[2];
sx q[2];
rz(-2.006861) q[2];
sx q[2];
rz(3.1046253) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0335281) q[1];
sx q[1];
rz(-2.2205097) q[1];
sx q[1];
rz(0.65945464) q[1];
x q[2];
rz(1.085161) q[3];
sx q[3];
rz(-1.3240459) q[3];
sx q[3];
rz(0.3890872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6923328) q[2];
sx q[2];
rz(-1.924694) q[2];
sx q[2];
rz(0.24162351) q[2];
rz(-0.062945098) q[3];
sx q[3];
rz(-1.2017622) q[3];
sx q[3];
rz(2.3285749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45074201) q[0];
sx q[0];
rz(-1.2797322) q[0];
sx q[0];
rz(-0.29020852) q[0];
rz(2.8065575) q[1];
sx q[1];
rz(-0.93828833) q[1];
sx q[1];
rz(0.32523528) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37595108) q[0];
sx q[0];
rz(-1.1055357) q[0];
sx q[0];
rz(-1.6164533) q[0];
rz(-pi) q[1];
rz(-1.8951072) q[2];
sx q[2];
rz(-0.70581573) q[2];
sx q[2];
rz(1.2845662) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.27694459) q[1];
sx q[1];
rz(-1.7200302) q[1];
sx q[1];
rz(1.4008888) q[1];
x q[2];
rz(0.70230647) q[3];
sx q[3];
rz(-0.78506535) q[3];
sx q[3];
rz(0.64083523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30514303) q[2];
sx q[2];
rz(-0.93905753) q[2];
sx q[2];
rz(-1.1391501) q[2];
rz(-2.3818453) q[3];
sx q[3];
rz(-0.097948827) q[3];
sx q[3];
rz(0.37794149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29001319) q[0];
sx q[0];
rz(-0.72999287) q[0];
sx q[0];
rz(2.3082025) q[0];
rz(-0.003412811) q[1];
sx q[1];
rz(-0.5178057) q[1];
sx q[1];
rz(1.980967) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7297555) q[0];
sx q[0];
rz(-1.5416391) q[0];
sx q[0];
rz(-1.2458318) q[0];
rz(2.4936475) q[2];
sx q[2];
rz(-1.750769) q[2];
sx q[2];
rz(-2.5729716) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6033481) q[1];
sx q[1];
rz(-1.5711492) q[1];
sx q[1];
rz(1.805549) q[1];
x q[2];
rz(-2.9594775) q[3];
sx q[3];
rz(-1.9312177) q[3];
sx q[3];
rz(-2.6153236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87355906) q[2];
sx q[2];
rz(-1.7256836) q[2];
sx q[2];
rz(2.9746919) q[2];
rz(1.9514826) q[3];
sx q[3];
rz(-0.24470617) q[3];
sx q[3];
rz(-2.0580097) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13389182) q[0];
sx q[0];
rz(-1.5950483) q[0];
sx q[0];
rz(0.85064763) q[0];
rz(0.8029241) q[1];
sx q[1];
rz(-1.9839169) q[1];
sx q[1];
rz(2.0096774) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9659289) q[0];
sx q[0];
rz(-0.53232876) q[0];
sx q[0];
rz(-1.2086297) q[0];
x q[1];
rz(1.8649624) q[2];
sx q[2];
rz(-1.7554211) q[2];
sx q[2];
rz(0.52900865) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2537576) q[1];
sx q[1];
rz(-1.3030392) q[1];
sx q[1];
rz(-2.5303809) q[1];
rz(-1.9798801) q[3];
sx q[3];
rz(-2.2581824) q[3];
sx q[3];
rz(-2.0076755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6159281) q[2];
sx q[2];
rz(-1.6388845) q[2];
sx q[2];
rz(-0.28044236) q[2];
rz(-0.40124714) q[3];
sx q[3];
rz(-0.29956996) q[3];
sx q[3];
rz(1.3569328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.930437) q[0];
sx q[0];
rz(-1.3933975) q[0];
sx q[0];
rz(2.9537971) q[0];
rz(0.64741778) q[1];
sx q[1];
rz(-0.96962601) q[1];
sx q[1];
rz(0.59026778) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3844073) q[0];
sx q[0];
rz(-1.2110855) q[0];
sx q[0];
rz(-3.1076184) q[0];
rz(-pi) q[1];
rz(-2.3529604) q[2];
sx q[2];
rz(-1.3071539) q[2];
sx q[2];
rz(-2.2396127) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62377159) q[1];
sx q[1];
rz(-1.54269) q[1];
sx q[1];
rz(0.080561056) q[1];
rz(0.26735683) q[3];
sx q[3];
rz(-0.68213576) q[3];
sx q[3];
rz(0.29058829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4831627) q[2];
sx q[2];
rz(-1.9715318) q[2];
sx q[2];
rz(-0.17407334) q[2];
rz(-0.46711323) q[3];
sx q[3];
rz(-0.73254782) q[3];
sx q[3];
rz(0.23021255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1751404) q[0];
sx q[0];
rz(-1.3850965) q[0];
sx q[0];
rz(2.0542282) q[0];
rz(0.1611791) q[1];
sx q[1];
rz(-1.9728856) q[1];
sx q[1];
rz(0.93343121) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1891834) q[0];
sx q[0];
rz(-1.2564474) q[0];
sx q[0];
rz(1.0745924) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1968139) q[2];
sx q[2];
rz(-1.6829964) q[2];
sx q[2];
rz(-2.2430132) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5111635) q[1];
sx q[1];
rz(-1.1463209) q[1];
sx q[1];
rz(-1.9349348) q[1];
rz(-2.3087738) q[3];
sx q[3];
rz(-0.76945451) q[3];
sx q[3];
rz(-0.11386816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33438385) q[2];
sx q[2];
rz(-1.4391359) q[2];
sx q[2];
rz(1.2558233) q[2];
rz(-3.1387591) q[3];
sx q[3];
rz(-2.5815559) q[3];
sx q[3];
rz(-0.66633666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46347076) q[0];
sx q[0];
rz(-3.1246298) q[0];
sx q[0];
rz(-0.28067881) q[0];
rz(0.30411389) q[1];
sx q[1];
rz(-0.76144832) q[1];
sx q[1];
rz(-1.4842518) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69347914) q[0];
sx q[0];
rz(-1.5511314) q[0];
sx q[0];
rz(-3.124191) q[0];
rz(0.78042729) q[2];
sx q[2];
rz(-1.4533416) q[2];
sx q[2];
rz(0.800363) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4485943) q[1];
sx q[1];
rz(-2.0830375) q[1];
sx q[1];
rz(-0.97657852) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1969332) q[3];
sx q[3];
rz(-1.6955293) q[3];
sx q[3];
rz(-2.9880808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41875276) q[2];
sx q[2];
rz(-0.80416983) q[2];
sx q[2];
rz(1.0214825) q[2];
rz(2.569765) q[3];
sx q[3];
rz(-0.56838667) q[3];
sx q[3];
rz(-2.2233326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9224213) q[0];
sx q[0];
rz(-2.2497441) q[0];
sx q[0];
rz(-0.047155596) q[0];
rz(-1.4157408) q[1];
sx q[1];
rz(-1.163131) q[1];
sx q[1];
rz(0.39173752) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3211818) q[0];
sx q[0];
rz(-1.5279378) q[0];
sx q[0];
rz(-2.6880797) q[0];
rz(-pi) q[1];
rz(1.0739296) q[2];
sx q[2];
rz(-2.7598682) q[2];
sx q[2];
rz(3.0860902) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.67591813) q[1];
sx q[1];
rz(-1.5603702) q[1];
sx q[1];
rz(-1.584076) q[1];
x q[2];
rz(2.6542386) q[3];
sx q[3];
rz(-0.4385786) q[3];
sx q[3];
rz(2.6233332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2262816) q[2];
sx q[2];
rz(-2.2690124) q[2];
sx q[2];
rz(2.937781) q[2];
rz(-2.5053744) q[3];
sx q[3];
rz(-1.3602942) q[3];
sx q[3];
rz(3.0740331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5930475) q[0];
sx q[0];
rz(-3.1199582) q[0];
sx q[0];
rz(-2.0245323) q[0];
rz(1.8695658) q[1];
sx q[1];
rz(-0.8465603) q[1];
sx q[1];
rz(2.5844432) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37763835) q[0];
sx q[0];
rz(-1.391811) q[0];
sx q[0];
rz(-1.1338439) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0821758) q[2];
sx q[2];
rz(-2.9520028) q[2];
sx q[2];
rz(-2.1185045) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0275201) q[1];
sx q[1];
rz(-1.2344242) q[1];
sx q[1];
rz(2.8739424) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8548202) q[3];
sx q[3];
rz(-1.0470812) q[3];
sx q[3];
rz(0.70543462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40994) q[2];
sx q[2];
rz(-0.2356379) q[2];
sx q[2];
rz(1.635599) q[2];
rz(-2.951494) q[3];
sx q[3];
rz(-0.59571576) q[3];
sx q[3];
rz(-3.0715517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53368038) q[0];
sx q[0];
rz(-0.30617014) q[0];
sx q[0];
rz(0.26468563) q[0];
rz(-2.7429122) q[1];
sx q[1];
rz(-0.52972263) q[1];
sx q[1];
rz(0.20464373) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5335352) q[0];
sx q[0];
rz(-1.7711955) q[0];
sx q[0];
rz(2.1222369) q[0];
x q[1];
rz(2.2024353) q[2];
sx q[2];
rz(-1.3654786) q[2];
sx q[2];
rz(3.1373295) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8890587) q[1];
sx q[1];
rz(-1.5100579) q[1];
sx q[1];
rz(-2.3142458) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4253623) q[3];
sx q[3];
rz(-2.0441291) q[3];
sx q[3];
rz(-2.4504599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9201422) q[2];
sx q[2];
rz(-1.7998989) q[2];
sx q[2];
rz(-1.7500925) q[2];
rz(-2.5642388) q[3];
sx q[3];
rz(-2.5632863) q[3];
sx q[3];
rz(0.81380832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6089384) q[0];
sx q[0];
rz(-1.895099) q[0];
sx q[0];
rz(1.1005703) q[0];
rz(-0.34250034) q[1];
sx q[1];
rz(-2.1771912) q[1];
sx q[1];
rz(2.049581) q[1];
rz(-1.2630812) q[2];
sx q[2];
rz(-1.7379496) q[2];
sx q[2];
rz(-1.3676277) q[2];
rz(-2.5534775) q[3];
sx q[3];
rz(-2.2089403) q[3];
sx q[3];
rz(-1.5939286) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
