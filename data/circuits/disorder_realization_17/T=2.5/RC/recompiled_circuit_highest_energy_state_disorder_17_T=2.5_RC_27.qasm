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
rz(1.0139581) q[0];
sx q[0];
rz(4.5944144) q[0];
sx q[0];
rz(8.8267427) q[0];
rz(1.1818089) q[1];
sx q[1];
rz(-0.66507566) q[1];
sx q[1];
rz(-2.9887078) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11422296) q[0];
sx q[0];
rz(-1.4589757) q[0];
sx q[0];
rz(-3.0614475) q[0];
rz(-pi) q[1];
rz(1.5122732) q[2];
sx q[2];
rz(-1.0835179) q[2];
sx q[2];
rz(-2.6807356) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8788556) q[1];
sx q[1];
rz(-1.3595743) q[1];
sx q[1];
rz(3.0131574) q[1];
rz(3.0769807) q[3];
sx q[3];
rz(-3.0955394) q[3];
sx q[3];
rz(-2.8833431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87024706) q[2];
sx q[2];
rz(-2.953244) q[2];
sx q[2];
rz(-1.9625473) q[2];
rz(2.7315268) q[3];
sx q[3];
rz(-1.054801) q[3];
sx q[3];
rz(2.2470391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.1882741) q[0];
sx q[0];
rz(-0.66903791) q[0];
sx q[0];
rz(2.8642995) q[0];
rz(-2.9412728) q[1];
sx q[1];
rz(-0.89626139) q[1];
sx q[1];
rz(-3.0576113) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5008299) q[0];
sx q[0];
rz(-2.6272054) q[0];
sx q[0];
rz(1.1577093) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7189156) q[2];
sx q[2];
rz(-2.2930055) q[2];
sx q[2];
rz(0.87004694) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9878432) q[1];
sx q[1];
rz(-1.1293518) q[1];
sx q[1];
rz(2.6325339) q[1];
x q[2];
rz(1.3091607) q[3];
sx q[3];
rz(-1.0519093) q[3];
sx q[3];
rz(-0.46552697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1818485) q[2];
sx q[2];
rz(-0.66755787) q[2];
sx q[2];
rz(0.7461156) q[2];
rz(-0.6212081) q[3];
sx q[3];
rz(-0.64380232) q[3];
sx q[3];
rz(1.752689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8212432) q[0];
sx q[0];
rz(-2.4140883) q[0];
sx q[0];
rz(-1.2028836) q[0];
rz(0.41162833) q[1];
sx q[1];
rz(-1.8609214) q[1];
sx q[1];
rz(-1.113755) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4899556) q[0];
sx q[0];
rz(-2.2764766) q[0];
sx q[0];
rz(-0.56244992) q[0];
rz(0.95543785) q[2];
sx q[2];
rz(-1.2005383) q[2];
sx q[2];
rz(0.40981612) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74677709) q[1];
sx q[1];
rz(-1.703843) q[1];
sx q[1];
rz(-0.89292553) q[1];
x q[2];
rz(2.2334072) q[3];
sx q[3];
rz(-0.81530967) q[3];
sx q[3];
rz(-0.96409982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0959629) q[2];
sx q[2];
rz(-0.7926422) q[2];
sx q[2];
rz(-1.1620129) q[2];
rz(0.80471936) q[3];
sx q[3];
rz(-1.2730803) q[3];
sx q[3];
rz(1.2812251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3941037) q[0];
sx q[0];
rz(-1.8590314) q[0];
sx q[0];
rz(-1.6266741) q[0];
rz(-0.50846848) q[1];
sx q[1];
rz(-0.6503121) q[1];
sx q[1];
rz(-1.633684) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8430458) q[0];
sx q[0];
rz(-0.15104476) q[0];
sx q[0];
rz(-1.5566948) q[0];
rz(-1.3541578) q[2];
sx q[2];
rz(-1.4329974) q[2];
sx q[2];
rz(2.0462284) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0993938) q[1];
sx q[1];
rz(-1.458199) q[1];
sx q[1];
rz(-1.8774752) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7135777) q[3];
sx q[3];
rz(-0.95967442) q[3];
sx q[3];
rz(-2.1111272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7278829) q[2];
sx q[2];
rz(-1.1943613) q[2];
sx q[2];
rz(-2.3940562) q[2];
rz(-1.2306635) q[3];
sx q[3];
rz(-2.3321798) q[3];
sx q[3];
rz(2.6443853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43664765) q[0];
sx q[0];
rz(-2.0132988) q[0];
sx q[0];
rz(-1.4592015) q[0];
rz(-2.7875426) q[1];
sx q[1];
rz(-0.67146462) q[1];
sx q[1];
rz(0.57463247) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.282772) q[0];
sx q[0];
rz(-3.0631815) q[0];
sx q[0];
rz(0.17228244) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28038763) q[2];
sx q[2];
rz(-0.78885733) q[2];
sx q[2];
rz(2.7352509) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3012817) q[1];
sx q[1];
rz(-1.5191829) q[1];
sx q[1];
rz(-0.2272931) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0784723) q[3];
sx q[3];
rz(-2.9083544) q[3];
sx q[3];
rz(2.4539029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.5057126) q[2];
sx q[2];
rz(-1.9680223) q[2];
sx q[2];
rz(2.9393348) q[2];
rz(-2.108719) q[3];
sx q[3];
rz(-2.5346916) q[3];
sx q[3];
rz(0.066298299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.1101087) q[0];
sx q[0];
rz(-2.7543572) q[0];
sx q[0];
rz(2.7701344) q[0];
rz(2.8455041) q[1];
sx q[1];
rz(-1.5541872) q[1];
sx q[1];
rz(-2.9782226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59054331) q[0];
sx q[0];
rz(-1.3871017) q[0];
sx q[0];
rz(-0.50409533) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4642015) q[2];
sx q[2];
rz(-2.4755726) q[2];
sx q[2];
rz(0.99413727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5886053) q[1];
sx q[1];
rz(-1.3901911) q[1];
sx q[1];
rz(2.6957979) q[1];
rz(-pi) q[2];
rz(1.7512239) q[3];
sx q[3];
rz(-0.70367763) q[3];
sx q[3];
rz(1.1168407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0614193) q[2];
sx q[2];
rz(-2.8714608) q[2];
sx q[2];
rz(-2.488625) q[2];
rz(2.669892) q[3];
sx q[3];
rz(-2.3931849) q[3];
sx q[3];
rz(0.72699839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93040526) q[0];
sx q[0];
rz(-1.0791595) q[0];
sx q[0];
rz(-2.1042673) q[0];
rz(-0.51517454) q[1];
sx q[1];
rz(-1.8637135) q[1];
sx q[1];
rz(-1.8151262) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4907287) q[0];
sx q[0];
rz(-2.0959637) q[0];
sx q[0];
rz(-2.690171) q[0];
rz(-pi) q[1];
rz(-1.5303361) q[2];
sx q[2];
rz(-1.9225112) q[2];
sx q[2];
rz(-0.88592285) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4660037) q[1];
sx q[1];
rz(-1.4741305) q[1];
sx q[1];
rz(-0.76670353) q[1];
rz(-0.063422261) q[3];
sx q[3];
rz(-1.2655971) q[3];
sx q[3];
rz(-2.2535107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4482164) q[2];
sx q[2];
rz(-0.46458149) q[2];
sx q[2];
rz(0.32619897) q[2];
rz(0.97801963) q[3];
sx q[3];
rz(-1.8947442) q[3];
sx q[3];
rz(3.0514362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4433032) q[0];
sx q[0];
rz(-2.8841618) q[0];
sx q[0];
rz(0.28939104) q[0];
rz(0.012271317) q[1];
sx q[1];
rz(-0.27996501) q[1];
sx q[1];
rz(1.4853005) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8809562) q[0];
sx q[0];
rz(-1.0614961) q[0];
sx q[0];
rz(-0.80672046) q[0];
rz(0.22356914) q[2];
sx q[2];
rz(-1.6235329) q[2];
sx q[2];
rz(1.2745672) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0721133) q[1];
sx q[1];
rz(-2.8918307) q[1];
sx q[1];
rz(0.62432557) q[1];
x q[2];
rz(2.0391123) q[3];
sx q[3];
rz(-1.2833929) q[3];
sx q[3];
rz(1.772955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5926008) q[2];
sx q[2];
rz(-1.5718549) q[2];
sx q[2];
rz(-2.9686417) q[2];
rz(-1.3744099) q[3];
sx q[3];
rz(-1.7632615) q[3];
sx q[3];
rz(-1.6636728) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7462815) q[0];
sx q[0];
rz(-0.44773856) q[0];
sx q[0];
rz(-0.81004274) q[0];
rz(2.7250302) q[1];
sx q[1];
rz(-0.77605334) q[1];
sx q[1];
rz(-0.3449482) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11703466) q[0];
sx q[0];
rz(-2.2753365) q[0];
sx q[0];
rz(-1.4092833) q[0];
rz(-pi) q[1];
rz(1.0428997) q[2];
sx q[2];
rz(-0.32020346) q[2];
sx q[2];
rz(-2.7803583) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7564063) q[1];
sx q[1];
rz(-2.029772) q[1];
sx q[1];
rz(-0.36151691) q[1];
rz(1.4515321) q[3];
sx q[3];
rz(-1.4550617) q[3];
sx q[3];
rz(1.7657464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.714146) q[2];
sx q[2];
rz(-0.88627187) q[2];
sx q[2];
rz(-0.084058849) q[2];
rz(0.90703026) q[3];
sx q[3];
rz(-1.7232938) q[3];
sx q[3];
rz(2.1840054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8091938) q[0];
sx q[0];
rz(-0.087662307) q[0];
sx q[0];
rz(0.3122538) q[0];
rz(-2.9088083) q[1];
sx q[1];
rz(-0.39118958) q[1];
sx q[1];
rz(1.8544474) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85111831) q[0];
sx q[0];
rz(-2.2016207) q[0];
sx q[0];
rz(-1.1680956) q[0];
rz(-pi) q[1];
rz(-1.5239117) q[2];
sx q[2];
rz(-0.93218902) q[2];
sx q[2];
rz(-2.435911) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5405266) q[1];
sx q[1];
rz(-0.52072111) q[1];
sx q[1];
rz(-2.0784318) q[1];
rz(-2.5165043) q[3];
sx q[3];
rz(-1.2660789) q[3];
sx q[3];
rz(-2.8947322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9707668) q[2];
sx q[2];
rz(-1.8584741) q[2];
sx q[2];
rz(-1.6440561) q[2];
rz(1.0981285) q[3];
sx q[3];
rz(-2.4427755) q[3];
sx q[3];
rz(2.6523377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1371786) q[0];
sx q[0];
rz(-1.351384) q[0];
sx q[0];
rz(2.2030892) q[0];
rz(1.8347523) q[1];
sx q[1];
rz(-0.88941457) q[1];
sx q[1];
rz(-0.79759146) q[1];
rz(-0.14456476) q[2];
sx q[2];
rz(-0.88695405) q[2];
sx q[2];
rz(0.6443413) q[2];
rz(-1.3908006) q[3];
sx q[3];
rz(-2.4672361) q[3];
sx q[3];
rz(1.5005021) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
