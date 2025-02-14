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
rz(-0.66752783) q[0];
sx q[0];
rz(-2.3025371) q[0];
sx q[0];
rz(2.7091205) q[0];
rz(-2.4617545) q[1];
sx q[1];
rz(-1.6814657) q[1];
sx q[1];
rz(0.7661933) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941152) q[0];
sx q[0];
rz(-1.9355398) q[0];
sx q[0];
rz(-0.095679521) q[0];
rz(-pi) q[1];
rz(-1.7518429) q[2];
sx q[2];
rz(-1.5403325) q[2];
sx q[2];
rz(-1.9150645) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4369504) q[1];
sx q[1];
rz(-2.5931103) q[1];
sx q[1];
rz(-1.8352527) q[1];
rz(-2.2965374) q[3];
sx q[3];
rz(-0.18394205) q[3];
sx q[3];
rz(-1.2181653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55521518) q[2];
sx q[2];
rz(-0.60474288) q[2];
sx q[2];
rz(-2.4797454) q[2];
rz(1.9072748) q[3];
sx q[3];
rz(-1.581278) q[3];
sx q[3];
rz(2.9890649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2374275) q[0];
sx q[0];
rz(-0.65218848) q[0];
sx q[0];
rz(2.7320614) q[0];
rz(-0.25853363) q[1];
sx q[1];
rz(-1.9347609) q[1];
sx q[1];
rz(-0.55005598) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4625699) q[0];
sx q[0];
rz(-1.4840811) q[0];
sx q[0];
rz(-2.0458512) q[0];
x q[1];
rz(-0.8895572) q[2];
sx q[2];
rz(-0.6830789) q[2];
sx q[2];
rz(-1.4448593) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.066795863) q[1];
sx q[1];
rz(-2.6535712) q[1];
sx q[1];
rz(2.0171432) q[1];
rz(-1.2004545) q[3];
sx q[3];
rz(-2.3562288) q[3];
sx q[3];
rz(1.0496275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6253769) q[2];
sx q[2];
rz(-1.7552525) q[2];
sx q[2];
rz(-0.59164444) q[2];
rz(-2.6243465) q[3];
sx q[3];
rz(-1.1040684) q[3];
sx q[3];
rz(-0.57466093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2194694) q[0];
sx q[0];
rz(-1.8485494) q[0];
sx q[0];
rz(-1.0009694) q[0];
rz(1.4283904) q[1];
sx q[1];
rz(-2.3826471) q[1];
sx q[1];
rz(-3.0975814) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6881268) q[0];
sx q[0];
rz(-1.1926873) q[0];
sx q[0];
rz(-0.72909683) q[0];
rz(0.83389177) q[2];
sx q[2];
rz(-1.5331037) q[2];
sx q[2];
rz(-0.99601638) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21979867) q[1];
sx q[1];
rz(-2.0646446) q[1];
sx q[1];
rz(-0.7864248) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8007565) q[3];
sx q[3];
rz(-2.965117) q[3];
sx q[3];
rz(-0.0214487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.22237805) q[2];
sx q[2];
rz(-1.5927477) q[2];
sx q[2];
rz(0.10981336) q[2];
rz(2.5114656) q[3];
sx q[3];
rz(-2.5870242) q[3];
sx q[3];
rz(-2.8540247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5146062) q[0];
sx q[0];
rz(-1.0543178) q[0];
sx q[0];
rz(-1.0940254) q[0];
rz(-2.6679954) q[1];
sx q[1];
rz(-2.4351599) q[1];
sx q[1];
rz(0.72023645) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7313774) q[0];
sx q[0];
rz(-2.5654554) q[0];
sx q[0];
rz(2.558484) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5237775) q[2];
sx q[2];
rz(-1.2628619) q[2];
sx q[2];
rz(-1.5297254) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2536972) q[1];
sx q[1];
rz(-1.869939) q[1];
sx q[1];
rz(2.3810785) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81497808) q[3];
sx q[3];
rz(-1.8061823) q[3];
sx q[3];
rz(-0.3334563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.56095162) q[2];
sx q[2];
rz(-1.4731984) q[2];
sx q[2];
rz(-2.9894323) q[2];
rz(1.2946607) q[3];
sx q[3];
rz(-1.7052238) q[3];
sx q[3];
rz(1.367182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0691567) q[0];
sx q[0];
rz(-2.5788588) q[0];
sx q[0];
rz(0.38134545) q[0];
rz(-2.5533679) q[1];
sx q[1];
rz(-1.4071848) q[1];
sx q[1];
rz(0.064195976) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9682024) q[0];
sx q[0];
rz(-1.9406835) q[0];
sx q[0];
rz(-1.3805519) q[0];
x q[1];
rz(2.5975294) q[2];
sx q[2];
rz(-0.69101101) q[2];
sx q[2];
rz(2.8008583) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.076157454) q[1];
sx q[1];
rz(-2.3373051) q[1];
sx q[1];
rz(1.1706513) q[1];
rz(-pi) q[2];
rz(1.3649356) q[3];
sx q[3];
rz(-1.6904217) q[3];
sx q[3];
rz(0.72577602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3885145) q[2];
sx q[2];
rz(-1.3103176) q[2];
sx q[2];
rz(-0.7927967) q[2];
rz(-3.0068398) q[3];
sx q[3];
rz(-2.5760791) q[3];
sx q[3];
rz(1.4904862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7029) q[0];
sx q[0];
rz(-1.5543992) q[0];
sx q[0];
rz(0.95293522) q[0];
rz(-1.8903525) q[1];
sx q[1];
rz(-1.6798991) q[1];
sx q[1];
rz(-0.23426637) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15794663) q[0];
sx q[0];
rz(-1.0704213) q[0];
sx q[0];
rz(-1.4051132) q[0];
x q[1];
rz(-1.9769922) q[2];
sx q[2];
rz(-1.147454) q[2];
sx q[2];
rz(0.51423954) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.837296) q[1];
sx q[1];
rz(-2.1526744) q[1];
sx q[1];
rz(-3.0051662) q[1];
x q[2];
rz(-1.686342) q[3];
sx q[3];
rz(-2.503814) q[3];
sx q[3];
rz(1.4664283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17144063) q[2];
sx q[2];
rz(-2.5415387) q[2];
sx q[2];
rz(2.1947412) q[2];
rz(1.6998477) q[3];
sx q[3];
rz(-1.8141247) q[3];
sx q[3];
rz(-3.0202151) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9117208) q[0];
sx q[0];
rz(-2.1389102) q[0];
sx q[0];
rz(0.3748931) q[0];
rz(-0.74294576) q[1];
sx q[1];
rz(-2.3187147) q[1];
sx q[1];
rz(-0.65623647) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0716055) q[0];
sx q[0];
rz(-2.2334369) q[0];
sx q[0];
rz(-1.9020984) q[0];
rz(-pi) q[1];
rz(2.1383275) q[2];
sx q[2];
rz(-2.1250688) q[2];
sx q[2];
rz(-0.90855937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69493964) q[1];
sx q[1];
rz(-0.57523721) q[1];
sx q[1];
rz(0.72845643) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5498059) q[3];
sx q[3];
rz(-1.7370708) q[3];
sx q[3];
rz(0.70032728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52834213) q[2];
sx q[2];
rz(-2.2684147) q[2];
sx q[2];
rz(0.64995107) q[2];
rz(2.2925099) q[3];
sx q[3];
rz(-0.63747469) q[3];
sx q[3];
rz(-1.953663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36986247) q[0];
sx q[0];
rz(-1.5700392) q[0];
sx q[0];
rz(0.318203) q[0];
rz(-2.8089306) q[1];
sx q[1];
rz(-0.55325621) q[1];
sx q[1];
rz(2.7645848) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4174667) q[0];
sx q[0];
rz(-1.2340349) q[0];
sx q[0];
rz(-0.048857852) q[0];
x q[1];
rz(-3.105806) q[2];
sx q[2];
rz(-2.0376959) q[2];
sx q[2];
rz(-0.70718471) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1313427) q[1];
sx q[1];
rz(-2.2325062) q[1];
sx q[1];
rz(-2.0493839) q[1];
x q[2];
rz(-0.11628233) q[3];
sx q[3];
rz(-2.0512037) q[3];
sx q[3];
rz(1.2670328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0074924) q[2];
sx q[2];
rz(-1.2780122) q[2];
sx q[2];
rz(-0.93097869) q[2];
rz(-1.603294) q[3];
sx q[3];
rz(-1.337991) q[3];
sx q[3];
rz(1.3163756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5285434) q[0];
sx q[0];
rz(-0.40754023) q[0];
sx q[0];
rz(2.7000309) q[0];
rz(-2.8764797) q[1];
sx q[1];
rz(-1.6997063) q[1];
sx q[1];
rz(-1.5107059) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1800507) q[0];
sx q[0];
rz(-1.7921711) q[0];
sx q[0];
rz(-0.44261296) q[0];
x q[1];
rz(2.1622545) q[2];
sx q[2];
rz(-1.854164) q[2];
sx q[2];
rz(1.3457553) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7889622) q[1];
sx q[1];
rz(-1.6704541) q[1];
sx q[1];
rz(-1.8046735) q[1];
rz(-pi) q[2];
rz(0.92075121) q[3];
sx q[3];
rz(-1.8691321) q[3];
sx q[3];
rz(-2.8118097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.6049217) q[2];
sx q[2];
rz(-1.8646381) q[2];
sx q[2];
rz(-2.6729909) q[2];
rz(-1.0118777) q[3];
sx q[3];
rz(-1.4230655) q[3];
sx q[3];
rz(1.0812149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.45379797) q[0];
sx q[0];
rz(-0.84306222) q[0];
sx q[0];
rz(1.4367562) q[0];
rz(0.94802481) q[1];
sx q[1];
rz(-1.3399622) q[1];
sx q[1];
rz(-2.8678133) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.439246) q[0];
sx q[0];
rz(-2.9986023) q[0];
sx q[0];
rz(-2.6197394) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1379399) q[2];
sx q[2];
rz(-0.64622245) q[2];
sx q[2];
rz(-0.60690676) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.791509) q[1];
sx q[1];
rz(-2.4285762) q[1];
sx q[1];
rz(2.3886695) q[1];
x q[2];
rz(-0.89694743) q[3];
sx q[3];
rz(-1.2123143) q[3];
sx q[3];
rz(-1.4523539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.580299) q[2];
sx q[2];
rz(-2.0686006) q[2];
sx q[2];
rz(2.6936074) q[2];
rz(-2.5042846) q[3];
sx q[3];
rz(-2.1154805) q[3];
sx q[3];
rz(0.89504939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7359903) q[0];
sx q[0];
rz(-1.6163106) q[0];
sx q[0];
rz(1.6769782) q[0];
rz(-2.6276656) q[1];
sx q[1];
rz(-0.75427873) q[1];
sx q[1];
rz(-2.8819081) q[1];
rz(2.9630255) q[2];
sx q[2];
rz(-2.6926276) q[2];
sx q[2];
rz(-1.2603029) q[2];
rz(-0.00058978422) q[3];
sx q[3];
rz(-1.9851709) q[3];
sx q[3];
rz(-2.6125534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
