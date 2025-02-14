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
rz(2.7338602) q[0];
sx q[0];
rz(-1.2091577) q[0];
sx q[0];
rz(-1.6785167) q[0];
rz(-2.8830124) q[1];
sx q[1];
rz(-1.2376031) q[1];
sx q[1];
rz(0.13374506) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27510724) q[0];
sx q[0];
rz(-1.5648399) q[0];
sx q[0];
rz(1.4136825) q[0];
rz(-pi) q[1];
rz(0.76579801) q[2];
sx q[2];
rz(-2.4054689) q[2];
sx q[2];
rz(0.36020261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7647296) q[1];
sx q[1];
rz(-1.2697744) q[1];
sx q[1];
rz(0.55913289) q[1];
x q[2];
rz(1.5695581) q[3];
sx q[3];
rz(-2.9887325) q[3];
sx q[3];
rz(-2.3521956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3818843) q[2];
sx q[2];
rz(-2.3298161) q[2];
sx q[2];
rz(-1.3105357) q[2];
rz(1.7740907) q[3];
sx q[3];
rz(-1.1279227) q[3];
sx q[3];
rz(-3.1254356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1644208) q[0];
sx q[0];
rz(-1.1559957) q[0];
sx q[0];
rz(-2.1261862) q[0];
rz(2.7482391) q[1];
sx q[1];
rz(-1.4721847) q[1];
sx q[1];
rz(-0.70158395) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6664826) q[0];
sx q[0];
rz(-2.283264) q[0];
sx q[0];
rz(-3.0300167) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.010741269) q[2];
sx q[2];
rz(-1.2678384) q[2];
sx q[2];
rz(-0.078127351) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2021705) q[1];
sx q[1];
rz(-1.8864453) q[1];
sx q[1];
rz(0.079074665) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2982992) q[3];
sx q[3];
rz(-0.95423428) q[3];
sx q[3];
rz(0.17184251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1141438) q[2];
sx q[2];
rz(-1.0237209) q[2];
sx q[2];
rz(2.0360428) q[2];
rz(2.8607199) q[3];
sx q[3];
rz(-2.5307145) q[3];
sx q[3];
rz(-2.9893379) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1059145) q[0];
sx q[0];
rz(-0.67258251) q[0];
sx q[0];
rz(2.0689082) q[0];
rz(-0.013280344) q[1];
sx q[1];
rz(-1.5417136) q[1];
sx q[1];
rz(-0.57659155) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6629366) q[0];
sx q[0];
rz(-2.8462704) q[0];
sx q[0];
rz(0.22462019) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0579975) q[2];
sx q[2];
rz(-1.7033249) q[2];
sx q[2];
rz(-2.8797163) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84725897) q[1];
sx q[1];
rz(-2.0461966) q[1];
sx q[1];
rz(-0.21902276) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9535096) q[3];
sx q[3];
rz(-0.59186223) q[3];
sx q[3];
rz(0.6955717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2836634) q[2];
sx q[2];
rz(-2.2495705) q[2];
sx q[2];
rz(-0.64209765) q[2];
rz(2.5765431) q[3];
sx q[3];
rz(-2.8665906) q[3];
sx q[3];
rz(-2.9716085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86244407) q[0];
sx q[0];
rz(-1.2553517) q[0];
sx q[0];
rz(-0.24236648) q[0];
rz(2.6077121) q[1];
sx q[1];
rz(-2.4568074) q[1];
sx q[1];
rz(-1.4885611) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3405995) q[0];
sx q[0];
rz(-1.4998319) q[0];
sx q[0];
rz(-0.59619716) q[0];
x q[1];
rz(-0.42885355) q[2];
sx q[2];
rz(-1.2786713) q[2];
sx q[2];
rz(-3.0512864) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3900534) q[1];
sx q[1];
rz(-1.281847) q[1];
sx q[1];
rz(-0.76559905) q[1];
x q[2];
rz(-1.6107992) q[3];
sx q[3];
rz(-1.8703331) q[3];
sx q[3];
rz(1.2114204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36294595) q[2];
sx q[2];
rz(-1.8774418) q[2];
sx q[2];
rz(1.470559) q[2];
rz(-0.10739022) q[3];
sx q[3];
rz(-1.8348179) q[3];
sx q[3];
rz(1.86514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11112467) q[0];
sx q[0];
rz(-2.680439) q[0];
sx q[0];
rz(-1.2029458) q[0];
rz(-0.89490926) q[1];
sx q[1];
rz(-1.9590961) q[1];
sx q[1];
rz(-2.9295909) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1231736) q[0];
sx q[0];
rz(-2.1668262) q[0];
sx q[0];
rz(-2.6072748) q[0];
x q[1];
rz(-2.6486829) q[2];
sx q[2];
rz(-2.5625336) q[2];
sx q[2];
rz(-2.2185203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4350076) q[1];
sx q[1];
rz(-0.75453484) q[1];
sx q[1];
rz(-0.88891397) q[1];
x q[2];
rz(-1.5552484) q[3];
sx q[3];
rz(-0.55335303) q[3];
sx q[3];
rz(1.6503594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8733069) q[2];
sx q[2];
rz(-2.3640859) q[2];
sx q[2];
rz(-3.0262465) q[2];
rz(0.98571944) q[3];
sx q[3];
rz(-1.7306354) q[3];
sx q[3];
rz(-2.7058097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2301521) q[0];
sx q[0];
rz(-1.6807012) q[0];
sx q[0];
rz(-2.4440515) q[0];
rz(-2.7279834) q[1];
sx q[1];
rz(-1.4613084) q[1];
sx q[1];
rz(0.70346171) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.441012) q[0];
sx q[0];
rz(-1.3003775) q[0];
sx q[0];
rz(-1.7187814) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.037961) q[2];
sx q[2];
rz(-0.67801266) q[2];
sx q[2];
rz(-2.968766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9338434) q[1];
sx q[1];
rz(-1.2022756) q[1];
sx q[1];
rz(0.50362419) q[1];
rz(-pi) q[2];
rz(2.892214) q[3];
sx q[3];
rz(-1.2802534) q[3];
sx q[3];
rz(-0.16502608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0387663) q[2];
sx q[2];
rz(-0.26831728) q[2];
sx q[2];
rz(-1.212567) q[2];
rz(0.28193998) q[3];
sx q[3];
rz(-1.599267) q[3];
sx q[3];
rz(2.6350002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3459699) q[0];
sx q[0];
rz(-0.75241929) q[0];
sx q[0];
rz(-0.84683007) q[0];
rz(-2.1961191) q[1];
sx q[1];
rz(-1.5739601) q[1];
sx q[1];
rz(-0.6792773) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1219471) q[0];
sx q[0];
rz(-1.6473734) q[0];
sx q[0];
rz(2.6640253) q[0];
rz(-0.34946673) q[2];
sx q[2];
rz(-0.65378705) q[2];
sx q[2];
rz(0.9330627) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0019466) q[1];
sx q[1];
rz(-1.0602385) q[1];
sx q[1];
rz(-2.9278048) q[1];
rz(-pi) q[2];
rz(-3.1222217) q[3];
sx q[3];
rz(-1.2865515) q[3];
sx q[3];
rz(-1.5721878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2519553) q[2];
sx q[2];
rz(-0.95459443) q[2];
sx q[2];
rz(1.6471222) q[2];
rz(0.55388081) q[3];
sx q[3];
rz(-1.5364105) q[3];
sx q[3];
rz(-1.7249829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81903356) q[0];
sx q[0];
rz(-2.366134) q[0];
sx q[0];
rz(-2.0490647) q[0];
rz(-2.20772) q[1];
sx q[1];
rz(-1.7364419) q[1];
sx q[1];
rz(-2.3563103) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0781495) q[0];
sx q[0];
rz(-2.5012996) q[0];
sx q[0];
rz(2.3904353) q[0];
rz(0.90964093) q[2];
sx q[2];
rz(-0.65654901) q[2];
sx q[2];
rz(-2.5971589) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.685647) q[1];
sx q[1];
rz(-0.37044493) q[1];
sx q[1];
rz(-0.73729314) q[1];
rz(2.0075032) q[3];
sx q[3];
rz(-0.99815449) q[3];
sx q[3];
rz(-2.9331895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2698722) q[2];
sx q[2];
rz(-1.4464658) q[2];
sx q[2];
rz(0.63228697) q[2];
rz(2.2271683) q[3];
sx q[3];
rz(-1.881003) q[3];
sx q[3];
rz(-2.9818592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81894994) q[0];
sx q[0];
rz(-0.85700789) q[0];
sx q[0];
rz(2.6981165) q[0];
rz(2.8387198) q[1];
sx q[1];
rz(-2.5587406) q[1];
sx q[1];
rz(-1.3076967) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0511868) q[0];
sx q[0];
rz(-1.7406347) q[0];
sx q[0];
rz(0.30533653) q[0];
x q[1];
rz(0.52438122) q[2];
sx q[2];
rz(-2.5705227) q[2];
sx q[2];
rz(-1.2829765) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.57932944) q[1];
sx q[1];
rz(-2.0585357) q[1];
sx q[1];
rz(0.30178782) q[1];
x q[2];
rz(1.0293209) q[3];
sx q[3];
rz(-2.6063237) q[3];
sx q[3];
rz(-0.9891555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.12094721) q[2];
sx q[2];
rz(-1.0409313) q[2];
sx q[2];
rz(0.61332235) q[2];
rz(0.81765085) q[3];
sx q[3];
rz(-0.48912564) q[3];
sx q[3];
rz(-2.8221655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15449512) q[0];
sx q[0];
rz(-1.3767865) q[0];
sx q[0];
rz(0.35926551) q[0];
rz(-2.477395) q[1];
sx q[1];
rz(-1.8134873) q[1];
sx q[1];
rz(-0.42116234) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0764685) q[0];
sx q[0];
rz(-1.2026005) q[0];
sx q[0];
rz(0.23875956) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0053275544) q[2];
sx q[2];
rz(-1.498328) q[2];
sx q[2];
rz(-0.42124149) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56257427) q[1];
sx q[1];
rz(-2.0236602) q[1];
sx q[1];
rz(-0.29154204) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99160414) q[3];
sx q[3];
rz(-2.0957855) q[3];
sx q[3];
rz(2.2158509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8757214) q[2];
sx q[2];
rz(-0.096780626) q[2];
sx q[2];
rz(-1.9270012) q[2];
rz(-2.752979) q[3];
sx q[3];
rz(-0.92258421) q[3];
sx q[3];
rz(1.0987561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6471967) q[0];
sx q[0];
rz(-1.5460486) q[0];
sx q[0];
rz(-1.8478951) q[0];
rz(-2.9248059) q[1];
sx q[1];
rz(-2.4516791) q[1];
sx q[1];
rz(-2.6286415) q[1];
rz(-2.4642261) q[2];
sx q[2];
rz(-1.2268885) q[2];
sx q[2];
rz(2.3033768) q[2];
rz(-2.6445845) q[3];
sx q[3];
rz(-2.6217032) q[3];
sx q[3];
rz(1.1006127) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
