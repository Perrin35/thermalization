OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.14804949) q[0];
sx q[0];
rz(-0.43819675) q[0];
sx q[0];
rz(10.155805) q[0];
rz(-2.464715) q[1];
sx q[1];
rz(-0.69898611) q[1];
sx q[1];
rz(-2.5875523) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0813538) q[0];
sx q[0];
rz(-2.0610272) q[0];
sx q[0];
rz(-2.1054563) q[0];
rz(-2.0238766) q[2];
sx q[2];
rz(-0.52342285) q[2];
sx q[2];
rz(-0.40016178) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4198534) q[1];
sx q[1];
rz(-1.5370779) q[1];
sx q[1];
rz(2.2583302) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6632313) q[3];
sx q[3];
rz(-2.195092) q[3];
sx q[3];
rz(2.9251298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.379091) q[2];
sx q[2];
rz(-1.5753626) q[2];
sx q[2];
rz(-0.16538922) q[2];
rz(-2.9108544) q[3];
sx q[3];
rz(-2.8985891) q[3];
sx q[3];
rz(0.86133426) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6623401) q[0];
sx q[0];
rz(-0.32821822) q[0];
sx q[0];
rz(1.3586556) q[0];
rz(-1.6183629) q[1];
sx q[1];
rz(-2.3871469) q[1];
sx q[1];
rz(2.3014136) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0531033) q[0];
sx q[0];
rz(-1.8376956) q[0];
sx q[0];
rz(1.0841814) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1112059) q[2];
sx q[2];
rz(-0.75558096) q[2];
sx q[2];
rz(0.43967512) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.45073918) q[1];
sx q[1];
rz(-2.2652341) q[1];
sx q[1];
rz(-1.9574375) q[1];
rz(-pi) q[2];
rz(2.3760017) q[3];
sx q[3];
rz(-1.4642707) q[3];
sx q[3];
rz(-1.2953143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3885865) q[2];
sx q[2];
rz(-1.2331839) q[2];
sx q[2];
rz(-0.93442717) q[2];
rz(1.2855444) q[3];
sx q[3];
rz(-2.3617187) q[3];
sx q[3];
rz(2.2540895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11338209) q[0];
sx q[0];
rz(-1.0209571) q[0];
sx q[0];
rz(1.5895948) q[0];
rz(1.7734843) q[1];
sx q[1];
rz(-1.2721456) q[1];
sx q[1];
rz(2.0210463) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6404868) q[0];
sx q[0];
rz(-1.2418134) q[0];
sx q[0];
rz(1.0803243) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0765614) q[2];
sx q[2];
rz(-0.49190822) q[2];
sx q[2];
rz(0.075752346) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.2255133) q[1];
sx q[1];
rz(-0.54359964) q[1];
sx q[1];
rz(-0.17982843) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2854101) q[3];
sx q[3];
rz(-2.719398) q[3];
sx q[3];
rz(2.5528922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0997194) q[2];
sx q[2];
rz(-0.21328829) q[2];
sx q[2];
rz(2.8495157) q[2];
rz(1.9000351) q[3];
sx q[3];
rz(-1.440666) q[3];
sx q[3];
rz(-0.45274538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9860155) q[0];
sx q[0];
rz(-0.41322511) q[0];
sx q[0];
rz(-2.7098932) q[0];
rz(2.9875634) q[1];
sx q[1];
rz(-1.5715716) q[1];
sx q[1];
rz(-2.0808751) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6018538) q[0];
sx q[0];
rz(-2.665747) q[0];
sx q[0];
rz(0.81625508) q[0];
x q[1];
rz(-2.9633386) q[2];
sx q[2];
rz(-1.3298103) q[2];
sx q[2];
rz(-0.18496938) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1460675) q[1];
sx q[1];
rz(-1.3473099) q[1];
sx q[1];
rz(1.2362012) q[1];
rz(-pi) q[2];
rz(-0.68935518) q[3];
sx q[3];
rz(-2.8971147) q[3];
sx q[3];
rz(2.1357128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3920307) q[2];
sx q[2];
rz(-1.531484) q[2];
sx q[2];
rz(0.57382601) q[2];
rz(-0.057403322) q[3];
sx q[3];
rz(-1.6618238) q[3];
sx q[3];
rz(-0.81006947) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5350128) q[0];
sx q[0];
rz(-0.25595328) q[0];
sx q[0];
rz(0.25273299) q[0];
rz(2.9622954) q[1];
sx q[1];
rz(-1.7959692) q[1];
sx q[1];
rz(-1.7326694) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5104467) q[0];
sx q[0];
rz(-0.24872336) q[0];
sx q[0];
rz(-1.876271) q[0];
rz(-1.9581355) q[2];
sx q[2];
rz(-1.9692067) q[2];
sx q[2];
rz(0.18926316) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9252569) q[1];
sx q[1];
rz(-2.3489423) q[1];
sx q[1];
rz(1.3986194) q[1];
rz(-0.5770251) q[3];
sx q[3];
rz(-0.60363704) q[3];
sx q[3];
rz(-1.9848422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1285642) q[2];
sx q[2];
rz(-1.3821673) q[2];
sx q[2];
rz(-1.9405091) q[2];
rz(-2.3342093) q[3];
sx q[3];
rz(-1.7816593) q[3];
sx q[3];
rz(2.0090296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0627237) q[0];
sx q[0];
rz(-1.6428592) q[0];
sx q[0];
rz(2.3003182) q[0];
rz(1.9367283) q[1];
sx q[1];
rz(-1.3179904) q[1];
sx q[1];
rz(-1.0296317) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1618274) q[0];
sx q[0];
rz(-2.0812391) q[0];
sx q[0];
rz(-3.1296041) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6252343) q[2];
sx q[2];
rz(-1.8022924) q[2];
sx q[2];
rz(-1.4524492) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8779422) q[1];
sx q[1];
rz(-1.7494838) q[1];
sx q[1];
rz(1.5787175) q[1];
x q[2];
rz(-1.383841) q[3];
sx q[3];
rz(-1.2724981) q[3];
sx q[3];
rz(1.8108071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.010217696) q[2];
sx q[2];
rz(-2.3848332) q[2];
sx q[2];
rz(2.6436464) q[2];
rz(0.52305269) q[3];
sx q[3];
rz(-1.7224576) q[3];
sx q[3];
rz(-0.12652346) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2987357) q[0];
sx q[0];
rz(-3.0099478) q[0];
sx q[0];
rz(-0.81197062) q[0];
rz(-0.89093351) q[1];
sx q[1];
rz(-1.9782601) q[1];
sx q[1];
rz(-0.99937159) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10103664) q[0];
sx q[0];
rz(-2.479739) q[0];
sx q[0];
rz(-0.7579114) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4677202) q[2];
sx q[2];
rz(-1.6666045) q[2];
sx q[2];
rz(0.16703781) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33617299) q[1];
sx q[1];
rz(-2.3567356) q[1];
sx q[1];
rz(-0.71342567) q[1];
rz(-2.2362367) q[3];
sx q[3];
rz(-1.3398583) q[3];
sx q[3];
rz(-0.11042102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.31819764) q[2];
sx q[2];
rz(-0.91788936) q[2];
sx q[2];
rz(1.8434175) q[2];
rz(0.038481742) q[3];
sx q[3];
rz(-1.1063856) q[3];
sx q[3];
rz(1.5255671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0865974) q[0];
sx q[0];
rz(-1.1308068) q[0];
sx q[0];
rz(2.8259592) q[0];
rz(1.2339969) q[1];
sx q[1];
rz(-0.71593586) q[1];
sx q[1];
rz(-1.2589781) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6427073) q[0];
sx q[0];
rz(-2.2872637) q[0];
sx q[0];
rz(1.1906719) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76541111) q[2];
sx q[2];
rz(-1.4452626) q[2];
sx q[2];
rz(1.5794992) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.087480672) q[1];
sx q[1];
rz(-1.9671983) q[1];
sx q[1];
rz(1.1328843) q[1];
rz(-pi) q[2];
rz(-3.0397814) q[3];
sx q[3];
rz(-1.5649127) q[3];
sx q[3];
rz(-2.7525097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9180318) q[2];
sx q[2];
rz(-2.3980902) q[2];
sx q[2];
rz(0.027916748) q[2];
rz(0.74808407) q[3];
sx q[3];
rz(-1.7223765) q[3];
sx q[3];
rz(2.9885651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41392031) q[0];
sx q[0];
rz(-0.99069178) q[0];
sx q[0];
rz(2.682611) q[0];
rz(2.0052295) q[1];
sx q[1];
rz(-1.7106979) q[1];
sx q[1];
rz(-2.9232025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55478375) q[0];
sx q[0];
rz(-3.112957) q[0];
sx q[0];
rz(0.61858004) q[0];
x q[1];
rz(1.893546) q[2];
sx q[2];
rz(-0.40229978) q[2];
sx q[2];
rz(-2.2108271) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4187516) q[1];
sx q[1];
rz(-1.7178078) q[1];
sx q[1];
rz(2.645869) q[1];
rz(-pi) q[2];
rz(1.9588542) q[3];
sx q[3];
rz(-1.6922631) q[3];
sx q[3];
rz(1.1282008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2961262) q[2];
sx q[2];
rz(-0.25456905) q[2];
sx q[2];
rz(-1.0477585) q[2];
rz(0.74538499) q[3];
sx q[3];
rz(-1.5785917) q[3];
sx q[3];
rz(1.4926225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1202241) q[0];
sx q[0];
rz(-0.57149082) q[0];
sx q[0];
rz(2.7987203) q[0];
rz(-2.951237) q[1];
sx q[1];
rz(-1.0089259) q[1];
sx q[1];
rz(-0.51317936) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1091217) q[0];
sx q[0];
rz(-2.4801804) q[0];
sx q[0];
rz(1.3107857) q[0];
rz(-pi) q[1];
rz(2.076976) q[2];
sx q[2];
rz(-1.2551885) q[2];
sx q[2];
rz(3.0402355) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4076297) q[1];
sx q[1];
rz(-1.9701951) q[1];
sx q[1];
rz(-0.49800378) q[1];
x q[2];
rz(-1.9916198) q[3];
sx q[3];
rz(-1.7003891) q[3];
sx q[3];
rz(-2.0146927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7036983) q[2];
sx q[2];
rz(-0.63592211) q[2];
sx q[2];
rz(-0.61100125) q[2];
rz(0.25887394) q[3];
sx q[3];
rz(-1.9301819) q[3];
sx q[3];
rz(1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7543058) q[0];
sx q[0];
rz(-1.5911234) q[0];
sx q[0];
rz(1.57244) q[0];
rz(0.98987956) q[1];
sx q[1];
rz(-1.2073333) q[1];
sx q[1];
rz(-2.5055199) q[1];
rz(0.0089916747) q[2];
sx q[2];
rz(-1.7594382) q[2];
sx q[2];
rz(0.3450763) q[2];
rz(1.5045139) q[3];
sx q[3];
rz(-2.2687329) q[3];
sx q[3];
rz(1.188059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
