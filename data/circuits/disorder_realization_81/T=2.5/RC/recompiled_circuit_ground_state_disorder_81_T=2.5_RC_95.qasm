OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17570198) q[0];
sx q[0];
rz(5.0134563) q[0];
sx q[0];
rz(11.391517) q[0];
rz(2.0570316) q[1];
sx q[1];
rz(4.6133981) q[1];
sx q[1];
rz(2.5785799) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63516894) q[0];
sx q[0];
rz(-0.20316589) q[0];
sx q[0];
rz(-0.26655339) q[0];
x q[1];
rz(-0.47663759) q[2];
sx q[2];
rz(-0.99299201) q[2];
sx q[2];
rz(3.1057242) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6463794) q[1];
sx q[1];
rz(-2.7732012) q[1];
sx q[1];
rz(0.45763335) q[1];
rz(-pi) q[2];
rz(2.7607542) q[3];
sx q[3];
rz(-1.9054753) q[3];
sx q[3];
rz(3.1370441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.53434831) q[2];
sx q[2];
rz(-1.4899985) q[2];
sx q[2];
rz(1.472817) q[2];
rz(1.8913174) q[3];
sx q[3];
rz(-1.6142802) q[3];
sx q[3];
rz(-0.97243029) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2431353) q[0];
sx q[0];
rz(-3.0861096) q[0];
sx q[0];
rz(2.6836416) q[0];
rz(-3.0398439) q[1];
sx q[1];
rz(-1.190217) q[1];
sx q[1];
rz(-2.8153458) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75263035) q[0];
sx q[0];
rz(-1.9174308) q[0];
sx q[0];
rz(-2.017157) q[0];
x q[1];
rz(-1.3497222) q[2];
sx q[2];
rz(-1.862163) q[2];
sx q[2];
rz(1.2779026) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.40405289) q[1];
sx q[1];
rz(-1.0657601) q[1];
sx q[1];
rz(2.2801599) q[1];
x q[2];
rz(-0.74298782) q[3];
sx q[3];
rz(-1.3338193) q[3];
sx q[3];
rz(-1.6176415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.132167) q[2];
sx q[2];
rz(-0.99401179) q[2];
sx q[2];
rz(1.2635292) q[2];
rz(-3.1215014) q[3];
sx q[3];
rz(-0.76954904) q[3];
sx q[3];
rz(-1.2300389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1175213) q[0];
sx q[0];
rz(-3.0776403) q[0];
sx q[0];
rz(0.57620311) q[0];
rz(1.4441215) q[1];
sx q[1];
rz(-1.2497808) q[1];
sx q[1];
rz(-1.261927) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44718868) q[0];
sx q[0];
rz(-1.2617636) q[0];
sx q[0];
rz(-0.033559994) q[0];
x q[1];
rz(-0.52719614) q[2];
sx q[2];
rz(-1.9771075) q[2];
sx q[2];
rz(-2.5715709) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.539725) q[1];
sx q[1];
rz(-0.66942184) q[1];
sx q[1];
rz(2.7407393) q[1];
rz(1.0141054) q[3];
sx q[3];
rz(-0.93645778) q[3];
sx q[3];
rz(-1.9012828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5835517) q[2];
sx q[2];
rz(-2.5463107) q[2];
sx q[2];
rz(-2.9884647) q[2];
rz(2.2554452) q[3];
sx q[3];
rz(-1.7821507) q[3];
sx q[3];
rz(-1.9396293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8837638) q[0];
sx q[0];
rz(-1.184329) q[0];
sx q[0];
rz(0.15876874) q[0];
rz(2.1067045) q[1];
sx q[1];
rz(-2.1885469) q[1];
sx q[1];
rz(0.75611702) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6203228) q[0];
sx q[0];
rz(-1.6286932) q[0];
sx q[0];
rz(1.6211364) q[0];
x q[1];
rz(1.4647843) q[2];
sx q[2];
rz(-2.8869625) q[2];
sx q[2];
rz(0.98843304) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0575916) q[1];
sx q[1];
rz(-2.822031) q[1];
sx q[1];
rz(2.4366385) q[1];
x q[2];
rz(1.3642374) q[3];
sx q[3];
rz(-1.8011923) q[3];
sx q[3];
rz(1.882706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5990344) q[2];
sx q[2];
rz(-0.90347806) q[2];
sx q[2];
rz(2.2912045) q[2];
rz(-1.4706069) q[3];
sx q[3];
rz(-1.8149523) q[3];
sx q[3];
rz(0.82408041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65664148) q[0];
sx q[0];
rz(-1.4018207) q[0];
sx q[0];
rz(-2.9630419) q[0];
rz(-1.8932331) q[1];
sx q[1];
rz(-1.5310042) q[1];
sx q[1];
rz(0.53370968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2399143) q[0];
sx q[0];
rz(-1.2216435) q[0];
sx q[0];
rz(1.3759744) q[0];
rz(0.99569624) q[2];
sx q[2];
rz(-0.786869) q[2];
sx q[2];
rz(2.5520419) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7099784) q[1];
sx q[1];
rz(-1.6233161) q[1];
sx q[1];
rz(2.2389328) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0664135) q[3];
sx q[3];
rz(-1.6565859) q[3];
sx q[3];
rz(2.498621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8278213) q[2];
sx q[2];
rz(-1.8905819) q[2];
sx q[2];
rz(2.0884464) q[2];
rz(0.18038067) q[3];
sx q[3];
rz(-0.80612055) q[3];
sx q[3];
rz(2.7789796) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5798222) q[0];
sx q[0];
rz(-2.1403911) q[0];
sx q[0];
rz(-1.9129821) q[0];
rz(-2.6749532) q[1];
sx q[1];
rz(-1.9231611) q[1];
sx q[1];
rz(-0.12900464) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0848964) q[0];
sx q[0];
rz(-1.7040623) q[0];
sx q[0];
rz(-0.10198089) q[0];
rz(0.84095793) q[2];
sx q[2];
rz(-1.691868) q[2];
sx q[2];
rz(0.88545495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1946757) q[1];
sx q[1];
rz(-1.2763775) q[1];
sx q[1];
rz(0.63398449) q[1];
rz(1.8390333) q[3];
sx q[3];
rz(-1.7058289) q[3];
sx q[3];
rz(2.6756833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1351607) q[2];
sx q[2];
rz(-2.5401523) q[2];
sx q[2];
rz(0.50430164) q[2];
rz(-2.048061) q[3];
sx q[3];
rz(-1.3239599) q[3];
sx q[3];
rz(-0.95660153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6495551) q[0];
sx q[0];
rz(-1.2832337) q[0];
sx q[0];
rz(-1.5851703) q[0];
rz(-0.57836142) q[1];
sx q[1];
rz(-1.2588986) q[1];
sx q[1];
rz(3.0392821) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033275398) q[0];
sx q[0];
rz(-0.96848291) q[0];
sx q[0];
rz(-2.6541416) q[0];
rz(-pi) q[1];
rz(-0.39686508) q[2];
sx q[2];
rz(-1.3584653) q[2];
sx q[2];
rz(0.020763606) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7150517) q[1];
sx q[1];
rz(-2.7965961) q[1];
sx q[1];
rz(2.662602) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3915406) q[3];
sx q[3];
rz(-2.1382209) q[3];
sx q[3];
rz(-0.22818434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1954701) q[2];
sx q[2];
rz(-1.2637694) q[2];
sx q[2];
rz(-0.61402399) q[2];
rz(0.92769235) q[3];
sx q[3];
rz(-1.5988908) q[3];
sx q[3];
rz(-0.6774925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(2.5880661) q[0];
sx q[0];
rz(-2.8480242) q[0];
sx q[0];
rz(-1.9422096) q[0];
rz(1.1792432) q[1];
sx q[1];
rz(-2.113138) q[1];
sx q[1];
rz(0.95338043) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48483585) q[0];
sx q[0];
rz(-2.0868882) q[0];
sx q[0];
rz(0.94658312) q[0];
x q[1];
rz(-3.0327263) q[2];
sx q[2];
rz(-1.035454) q[2];
sx q[2];
rz(-1.9809198) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.17705616) q[1];
sx q[1];
rz(-2.3971618) q[1];
sx q[1];
rz(-1.4076648) q[1];
rz(-pi) q[2];
rz(2.6598831) q[3];
sx q[3];
rz(-1.5007956) q[3];
sx q[3];
rz(2.9298129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.74932468) q[2];
sx q[2];
rz(-3.0159123) q[2];
sx q[2];
rz(-2.9607062) q[2];
rz(-2.0285897) q[3];
sx q[3];
rz(-1.0180611) q[3];
sx q[3];
rz(-0.38016144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.0522633) q[0];
sx q[0];
rz(-0.90173975) q[0];
sx q[0];
rz(2.3774636) q[0];
rz(1.0230505) q[1];
sx q[1];
rz(-1.5307531) q[1];
sx q[1];
rz(-1.4952362) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4039072) q[0];
sx q[0];
rz(-1.7244743) q[0];
sx q[0];
rz(-2.840848) q[0];
rz(-pi) q[1];
rz(2.1713397) q[2];
sx q[2];
rz(-1.2926599) q[2];
sx q[2];
rz(0.32626611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5896411) q[1];
sx q[1];
rz(-1.9560588) q[1];
sx q[1];
rz(-1.4937708) q[1];
rz(2.6999254) q[3];
sx q[3];
rz(-1.212647) q[3];
sx q[3];
rz(1.6902069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19227795) q[2];
sx q[2];
rz(-2.1239026) q[2];
sx q[2];
rz(2.6940572) q[2];
rz(-3.0053075) q[3];
sx q[3];
rz(-2.648573) q[3];
sx q[3];
rz(-0.59806699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78551453) q[0];
sx q[0];
rz(-1.0250174) q[0];
sx q[0];
rz(-2.6688975) q[0];
rz(-2.7418432) q[1];
sx q[1];
rz(-2.4374297) q[1];
sx q[1];
rz(-0.8185111) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8560573) q[0];
sx q[0];
rz(-1.7489079) q[0];
sx q[0];
rz(0.31383236) q[0];
rz(-pi) q[1];
rz(-1.1498951) q[2];
sx q[2];
rz(-2.2031261) q[2];
sx q[2];
rz(2.9078751) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.0014035066) q[1];
sx q[1];
rz(-1.2480987) q[1];
sx q[1];
rz(3.138209) q[1];
x q[2];
rz(2.6916041) q[3];
sx q[3];
rz(-2.7594824) q[3];
sx q[3];
rz(-1.9616322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58273903) q[2];
sx q[2];
rz(-2.0595198) q[2];
sx q[2];
rz(-1.2782798) q[2];
rz(1.9418955) q[3];
sx q[3];
rz(-1.4408828) q[3];
sx q[3];
rz(-2.8619158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3661135) q[0];
sx q[0];
rz(-1.7178602) q[0];
sx q[0];
rz(-1.3445509) q[0];
rz(1.7300425) q[1];
sx q[1];
rz(-0.81137864) q[1];
sx q[1];
rz(2.484533) q[1];
rz(3.1362202) q[2];
sx q[2];
rz(-3.0006144) q[2];
sx q[2];
rz(-0.49868546) q[2];
rz(1.1468519) q[3];
sx q[3];
rz(-0.70953555) q[3];
sx q[3];
rz(-2.4279688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
