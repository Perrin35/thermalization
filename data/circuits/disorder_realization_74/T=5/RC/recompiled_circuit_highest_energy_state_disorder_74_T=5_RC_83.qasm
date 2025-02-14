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
rz(-0.44218818) q[0];
sx q[0];
rz(-1.1531885) q[0];
sx q[0];
rz(-0.12864223) q[0];
rz(1.6755942) q[1];
sx q[1];
rz(-2.0425551) q[1];
sx q[1];
rz(1.0817945) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32444123) q[0];
sx q[0];
rz(-1.5400346) q[0];
sx q[0];
rz(-0.064662393) q[0];
x q[1];
rz(-2.2954582) q[2];
sx q[2];
rz(-0.48389176) q[2];
sx q[2];
rz(-1.0852551) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1638559) q[1];
sx q[1];
rz(-2.9675936) q[1];
sx q[1];
rz(-1.15088) q[1];
rz(1.2461189) q[3];
sx q[3];
rz(-0.28936181) q[3];
sx q[3];
rz(-2.112767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.094540207) q[2];
sx q[2];
rz(-2.4461942) q[2];
sx q[2];
rz(0.73235861) q[2];
rz(0.098946027) q[3];
sx q[3];
rz(-2.1061335) q[3];
sx q[3];
rz(-1.8100479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.052213) q[0];
sx q[0];
rz(-0.78240028) q[0];
sx q[0];
rz(-3.1001477) q[0];
rz(0.6913569) q[1];
sx q[1];
rz(-1.3203011) q[1];
sx q[1];
rz(2.2533805) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23032886) q[0];
sx q[0];
rz(-2.3155876) q[0];
sx q[0];
rz(-2.2780096) q[0];
rz(2.0574548) q[2];
sx q[2];
rz(-2.1900822) q[2];
sx q[2];
rz(-2.3145683) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8543635) q[1];
sx q[1];
rz(-1.3280265) q[1];
sx q[1];
rz(1.6153264) q[1];
rz(2.668374) q[3];
sx q[3];
rz(-2.5788973) q[3];
sx q[3];
rz(2.0317047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9899675) q[2];
sx q[2];
rz(-1.7513559) q[2];
sx q[2];
rz(-1.7612696) q[2];
rz(3.0917998) q[3];
sx q[3];
rz(-2.4536665) q[3];
sx q[3];
rz(-0.5602347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3675073) q[0];
sx q[0];
rz(-0.15223509) q[0];
sx q[0];
rz(1.2138858) q[0];
rz(1.3146575) q[1];
sx q[1];
rz(-1.0379125) q[1];
sx q[1];
rz(-1.5705869) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1477858) q[0];
sx q[0];
rz(-1.3249614) q[0];
sx q[0];
rz(-2.7509806) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7271829) q[2];
sx q[2];
rz(-1.1646574) q[2];
sx q[2];
rz(1.0473503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27109001) q[1];
sx q[1];
rz(-1.5368286) q[1];
sx q[1];
rz(-1.7351331) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96974428) q[3];
sx q[3];
rz(-1.0625524) q[3];
sx q[3];
rz(0.54304063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7870002) q[2];
sx q[2];
rz(-2.3921693) q[2];
sx q[2];
rz(-0.3053537) q[2];
rz(-2.2954156) q[3];
sx q[3];
rz(-1.927522) q[3];
sx q[3];
rz(-0.86887104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884035) q[0];
sx q[0];
rz(-2.1222293) q[0];
sx q[0];
rz(2.1715721) q[0];
rz(0.24523973) q[1];
sx q[1];
rz(-1.0673362) q[1];
sx q[1];
rz(1.6397569) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1872699) q[0];
sx q[0];
rz(-1.6980972) q[0];
sx q[0];
rz(-2.9856647) q[0];
x q[1];
rz(-0.63070339) q[2];
sx q[2];
rz(-2.8526488) q[2];
sx q[2];
rz(0.17221355) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40250999) q[1];
sx q[1];
rz(-1.7231047) q[1];
sx q[1];
rz(-2.9500752) q[1];
rz(-pi) q[2];
rz(2.4145441) q[3];
sx q[3];
rz(-1.8173976) q[3];
sx q[3];
rz(-1.5545157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.60845) q[2];
sx q[2];
rz(-1.5641944) q[2];
sx q[2];
rz(-1.8428141) q[2];
rz(-2.6070144) q[3];
sx q[3];
rz(-0.60408533) q[3];
sx q[3];
rz(-0.66876137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0534441) q[0];
sx q[0];
rz(-1.3548594) q[0];
sx q[0];
rz(0.97766367) q[0];
rz(2.4560302) q[1];
sx q[1];
rz(-2.3473163) q[1];
sx q[1];
rz(-0.96644863) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053542972) q[0];
sx q[0];
rz(-0.85659638) q[0];
sx q[0];
rz(2.5122253) q[0];
rz(-pi) q[1];
rz(0.62018779) q[2];
sx q[2];
rz(-0.66588565) q[2];
sx q[2];
rz(2.263139) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9453166) q[1];
sx q[1];
rz(-1.4320201) q[1];
sx q[1];
rz(0.86244418) q[1];
rz(2.590476) q[3];
sx q[3];
rz(-1.8752021) q[3];
sx q[3];
rz(-1.1486883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.55296772) q[2];
sx q[2];
rz(-1.3543465) q[2];
sx q[2];
rz(-2.8589613) q[2];
rz(0.96453729) q[3];
sx q[3];
rz(-1.5610361) q[3];
sx q[3];
rz(-0.36981043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.033217) q[0];
sx q[0];
rz(-2.7300457) q[0];
sx q[0];
rz(1.0785528) q[0];
rz(0.82303965) q[1];
sx q[1];
rz(-2.0185202) q[1];
sx q[1];
rz(-2.3413234) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5259994) q[0];
sx q[0];
rz(-1.5416) q[0];
sx q[0];
rz(0.008511624) q[0];
rz(0.79485017) q[2];
sx q[2];
rz(-2.7140116) q[2];
sx q[2];
rz(-2.7045369) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9046968) q[1];
sx q[1];
rz(-2.176732) q[1];
sx q[1];
rz(-0.69463457) q[1];
rz(-2.0144281) q[3];
sx q[3];
rz(-2.1882957) q[3];
sx q[3];
rz(-1.3408183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0018953) q[2];
sx q[2];
rz(-1.8997833) q[2];
sx q[2];
rz(-0.099420698) q[2];
rz(2.5771778) q[3];
sx q[3];
rz(-1.5494917) q[3];
sx q[3];
rz(0.91993371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99892202) q[0];
sx q[0];
rz(-0.26695928) q[0];
sx q[0];
rz(-0.024918407) q[0];
rz(1.092356) q[1];
sx q[1];
rz(-1.151029) q[1];
sx q[1];
rz(2.5999462) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6059581) q[0];
sx q[0];
rz(-2.1660457) q[0];
sx q[0];
rz(0.43979494) q[0];
rz(-1.0007896) q[2];
sx q[2];
rz(-1.2932245) q[2];
sx q[2];
rz(-2.9961207) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4389651) q[1];
sx q[1];
rz(-2.781233) q[1];
sx q[1];
rz(1.9907336) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3176407) q[3];
sx q[3];
rz(-2.5196155) q[3];
sx q[3];
rz(-0.013381026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8267374) q[2];
sx q[2];
rz(-1.22437) q[2];
sx q[2];
rz(-1.1401736) q[2];
rz(1.4402116) q[3];
sx q[3];
rz(-2.7482016) q[3];
sx q[3];
rz(2.9668729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88406554) q[0];
sx q[0];
rz(-1.9767569) q[0];
sx q[0];
rz(0.87673941) q[0];
rz(0.17403099) q[1];
sx q[1];
rz(-1.0935676) q[1];
sx q[1];
rz(1.7736951) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4583261) q[0];
sx q[0];
rz(-1.2948827) q[0];
sx q[0];
rz(2.1349195) q[0];
rz(-pi) q[1];
rz(2.9139745) q[2];
sx q[2];
rz(-1.926486) q[2];
sx q[2];
rz(1.571589) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0292082) q[1];
sx q[1];
rz(-1.9416326) q[1];
sx q[1];
rz(0.98354323) q[1];
rz(-0.075966751) q[3];
sx q[3];
rz(-1.6923762) q[3];
sx q[3];
rz(2.0437123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4725264) q[2];
sx q[2];
rz(-1.8146699) q[2];
sx q[2];
rz(0.18829045) q[2];
rz(-1.7752198) q[3];
sx q[3];
rz(-1.7732737) q[3];
sx q[3];
rz(-1.9395456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5240204) q[0];
sx q[0];
rz(-3.0078648) q[0];
sx q[0];
rz(-0.66620859) q[0];
rz(-2.0062402) q[1];
sx q[1];
rz(-2.1609781) q[1];
sx q[1];
rz(-1.1599783) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0157601) q[0];
sx q[0];
rz(-2.1938132) q[0];
sx q[0];
rz(-2.4218904) q[0];
x q[1];
rz(2.9097144) q[2];
sx q[2];
rz(-0.7674976) q[2];
sx q[2];
rz(-0.19814834) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.682232) q[1];
sx q[1];
rz(-2.5256093) q[1];
sx q[1];
rz(-1.5840696) q[1];
rz(-pi) q[2];
rz(-2.1362392) q[3];
sx q[3];
rz(-1.8400888) q[3];
sx q[3];
rz(2.5776742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3975415) q[2];
sx q[2];
rz(-0.56637374) q[2];
sx q[2];
rz(-0.19691021) q[2];
rz(2.2714553) q[3];
sx q[3];
rz(-1.0715485) q[3];
sx q[3];
rz(-1.8740694) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63721913) q[0];
sx q[0];
rz(-1.0727896) q[0];
sx q[0];
rz(-3.1043501) q[0];
rz(1.0875018) q[1];
sx q[1];
rz(-1.0481513) q[1];
sx q[1];
rz(-2.9748999) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0282306) q[0];
sx q[0];
rz(-0.89269243) q[0];
sx q[0];
rz(1.9166758) q[0];
rz(-pi) q[1];
rz(2.6564381) q[2];
sx q[2];
rz(-1.8388766) q[2];
sx q[2];
rz(1.9659276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8240451) q[1];
sx q[1];
rz(-2.8338441) q[1];
sx q[1];
rz(-2.8466562) q[1];
rz(-0.44346614) q[3];
sx q[3];
rz(-0.16032585) q[3];
sx q[3];
rz(-0.043930862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7249666) q[2];
sx q[2];
rz(-0.61369696) q[2];
sx q[2];
rz(1.7566768) q[2];
rz(0.40063217) q[3];
sx q[3];
rz(-1.2847565) q[3];
sx q[3];
rz(-1.0111151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2597802) q[0];
sx q[0];
rz(-2.0406944) q[0];
sx q[0];
rz(2.5115321) q[0];
rz(0.13314816) q[1];
sx q[1];
rz(-1.9825736) q[1];
sx q[1];
rz(2.0864743) q[1];
rz(-1.7693949) q[2];
sx q[2];
rz(-2.9469941) q[2];
sx q[2];
rz(-0.17175248) q[2];
rz(1.6937485) q[3];
sx q[3];
rz(-2.3364441) q[3];
sx q[3];
rz(0.91329109) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
