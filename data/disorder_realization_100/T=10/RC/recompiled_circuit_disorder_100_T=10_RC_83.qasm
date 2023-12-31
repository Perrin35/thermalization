OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(2.7913845) q[0];
sx q[0];
rz(12.933001) q[0];
rz(0.86759138) q[1];
sx q[1];
rz(-2.4974513) q[1];
sx q[1];
rz(-1.6860513) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.769387) q[0];
sx q[0];
rz(-1.1507478) q[0];
sx q[0];
rz(0.37680349) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26379649) q[2];
sx q[2];
rz(-2.2552239) q[2];
sx q[2];
rz(-0.85927187) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.75016025) q[1];
sx q[1];
rz(-0.62392985) q[1];
sx q[1];
rz(2.0425914) q[1];
x q[2];
rz(3.0936712) q[3];
sx q[3];
rz(-1.5600292) q[3];
sx q[3];
rz(-0.79706942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.1047487) q[2];
sx q[2];
rz(-1.3354744) q[2];
sx q[2];
rz(0.39164266) q[2];
rz(-0.13970217) q[3];
sx q[3];
rz(-0.67290664) q[3];
sx q[3];
rz(-0.02123775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7249107) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(-1.2765983) q[0];
rz(2.2712767) q[1];
sx q[1];
rz(-1.5753997) q[1];
sx q[1];
rz(1.93719) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6635839) q[0];
sx q[0];
rz(-1.1482129) q[0];
sx q[0];
rz(-0.43310662) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3538829) q[2];
sx q[2];
rz(-0.7173983) q[2];
sx q[2];
rz(-2.4153828) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7960297) q[1];
sx q[1];
rz(-1.2064762) q[1];
sx q[1];
rz(-1.6771392) q[1];
x q[2];
rz(2.2435917) q[3];
sx q[3];
rz(-1.3452531) q[3];
sx q[3];
rz(-2.0190092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5045972) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(1.2191999) q[2];
rz(0.35456625) q[3];
sx q[3];
rz(-0.48186007) q[3];
sx q[3];
rz(-1.5725117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59042674) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(2.5235126) q[0];
rz(0.016050054) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(-1.9504257) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1026099) q[0];
sx q[0];
rz(-1.8840186) q[0];
sx q[0];
rz(-0.15977504) q[0];
x q[1];
rz(-2.8580655) q[2];
sx q[2];
rz(-1.8054188) q[2];
sx q[2];
rz(1.7262176) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.788113) q[1];
sx q[1];
rz(-2.8444926) q[1];
sx q[1];
rz(-2.3658845) q[1];
x q[2];
rz(-0.61471625) q[3];
sx q[3];
rz(-0.63496548) q[3];
sx q[3];
rz(0.8599417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96737343) q[2];
sx q[2];
rz(-2.1790395) q[2];
sx q[2];
rz(-1.013914) q[2];
rz(1.7845456) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(-2.1092265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8106666) q[0];
sx q[0];
rz(-1.7174915) q[0];
sx q[0];
rz(0.20430918) q[0];
rz(1.7640242) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(-2.2185982) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5343691) q[0];
sx q[0];
rz(-1.9019706) q[0];
sx q[0];
rz(1.2481199) q[0];
rz(-pi) q[1];
x q[1];
rz(0.043945233) q[2];
sx q[2];
rz(-1.5825795) q[2];
sx q[2];
rz(-2.1859283) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1624565) q[1];
sx q[1];
rz(-0.41185954) q[1];
sx q[1];
rz(1.8331752) q[1];
x q[2];
rz(0.8664341) q[3];
sx q[3];
rz(-1.5537062) q[3];
sx q[3];
rz(0.41364241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4541645) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(0.50951177) q[2];
rz(-0.64309684) q[3];
sx q[3];
rz(-1.0984848) q[3];
sx q[3];
rz(-2.174214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.61280695) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(-0.91941419) q[0];
rz(-0.98584229) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(-1.7808419) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0589941) q[0];
sx q[0];
rz(-1.7599802) q[0];
sx q[0];
rz(-0.22342213) q[0];
rz(-0.947457) q[2];
sx q[2];
rz(-0.91295487) q[2];
sx q[2];
rz(0.62589494) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4701925) q[1];
sx q[1];
rz(-2.6365286) q[1];
sx q[1];
rz(0.43882521) q[1];
rz(-pi) q[2];
rz(-0.97134437) q[3];
sx q[3];
rz(-1.506862) q[3];
sx q[3];
rz(-1.5980699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78648606) q[2];
sx q[2];
rz(-1.3856709) q[2];
sx q[2];
rz(-0.11486593) q[2];
rz(-0.45587513) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(1.8615287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8961261) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(1.2448357) q[0];
rz(2.4339829) q[1];
sx q[1];
rz(-1.6376303) q[1];
sx q[1];
rz(0.27522603) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9305785) q[0];
sx q[0];
rz(-2.0654581) q[0];
sx q[0];
rz(-2.7633694) q[0];
rz(0.92991021) q[2];
sx q[2];
rz(-1.9800131) q[2];
sx q[2];
rz(0.90604679) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3616398) q[1];
sx q[1];
rz(-2.0556903) q[1];
sx q[1];
rz(1.1760902) q[1];
x q[2];
rz(3.0589468) q[3];
sx q[3];
rz(-2.189889) q[3];
sx q[3];
rz(1.9649399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.39367166) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(-2.6829524) q[2];
rz(2.4300872) q[3];
sx q[3];
rz(-2.4578874) q[3];
sx q[3];
rz(0.59035629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8608619) q[0];
sx q[0];
rz(-1.1870528) q[0];
sx q[0];
rz(1.77805) q[0];
rz(1.7215615) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(-0.71969676) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5170383) q[0];
sx q[0];
rz(-1.8554243) q[0];
sx q[0];
rz(0.65529234) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4160956) q[2];
sx q[2];
rz(-2.2273387) q[2];
sx q[2];
rz(-0.072349116) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3230348) q[1];
sx q[1];
rz(-2.5948988) q[1];
sx q[1];
rz(0.085303765) q[1];
rz(-pi) q[2];
rz(0.32960524) q[3];
sx q[3];
rz(-2.4334987) q[3];
sx q[3];
rz(3.0081188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.65934962) q[2];
sx q[2];
rz(-1.6094001) q[2];
sx q[2];
rz(2.4534295) q[2];
rz(-3.0996389) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(0.28013128) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91350895) q[0];
sx q[0];
rz(-1.1627473) q[0];
sx q[0];
rz(2.9550609) q[0];
rz(-0.60449156) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(-1.6465181) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1218425) q[0];
sx q[0];
rz(-1.9103721) q[0];
sx q[0];
rz(1.7975438) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69627701) q[2];
sx q[2];
rz(-1.7211203) q[2];
sx q[2];
rz(2.6339649) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1868134) q[1];
sx q[1];
rz(-0.28007945) q[1];
sx q[1];
rz(-2.7571452) q[1];
x q[2];
rz(2.2349615) q[3];
sx q[3];
rz(-2.2615221) q[3];
sx q[3];
rz(-0.074113473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.24215332) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(-3.126826) q[2];
rz(-1.8642558) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(-1.4130672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(2.9799141) q[0];
sx q[0];
rz(-2.4665687) q[0];
sx q[0];
rz(2.263608) q[0];
rz(-2.9453078) q[1];
sx q[1];
rz(-1.9613962) q[1];
sx q[1];
rz(1.3605798) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.774705) q[0];
sx q[0];
rz(-2.5686712) q[0];
sx q[0];
rz(1.6402871) q[0];
rz(-3.0948823) q[2];
sx q[2];
rz(-2.197406) q[2];
sx q[2];
rz(-1.4057297) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1189551) q[1];
sx q[1];
rz(-1.2345018) q[1];
sx q[1];
rz(-2.789546) q[1];
rz(-pi) q[2];
rz(1.9079886) q[3];
sx q[3];
rz(-0.4412776) q[3];
sx q[3];
rz(-2.6233167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62375162) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(2.2488135) q[2];
rz(3.1023074) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(0.75004309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84353012) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(0.046534006) q[0];
rz(0.90905601) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(-2.4216901) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35346183) q[0];
sx q[0];
rz(-2.8849368) q[0];
sx q[0];
rz(1.2528332) q[0];
rz(0.30883046) q[2];
sx q[2];
rz(-1.575982) q[2];
sx q[2];
rz(1.3877102) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0278922) q[1];
sx q[1];
rz(-0.24735951) q[1];
sx q[1];
rz(-1.5075831) q[1];
rz(0.87749691) q[3];
sx q[3];
rz(-0.96713669) q[3];
sx q[3];
rz(-2.7503777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7211192) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(-0.026467888) q[2];
rz(1.3039533) q[3];
sx q[3];
rz(-0.81726685) q[3];
sx q[3];
rz(-1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(1.3532886) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(2.9251255) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(-1.2555867) q[2];
sx q[2];
rz(-2.253058) q[2];
sx q[2];
rz(-0.81894973) q[2];
rz(-3.1090267) q[3];
sx q[3];
rz(-2.1571772) q[3];
sx q[3];
rz(-1.0909506) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
