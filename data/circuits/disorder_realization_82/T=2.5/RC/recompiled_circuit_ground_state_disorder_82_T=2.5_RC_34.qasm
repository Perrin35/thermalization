OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7656443) q[0];
sx q[0];
rz(-0.41416895) q[0];
sx q[0];
rz(2.3400657) q[0];
rz(3.1385359) q[1];
sx q[1];
rz(-2.2771775) q[1];
sx q[1];
rz(-3.0473696) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021359062) q[0];
sx q[0];
rz(-0.84478837) q[0];
sx q[0];
rz(1.057522) q[0];
rz(-pi) q[1];
rz(-3.1108833) q[2];
sx q[2];
rz(-1.2982663) q[2];
sx q[2];
rz(-2.3365793) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4068216) q[1];
sx q[1];
rz(-2.0367921) q[1];
sx q[1];
rz(-2.8193406) q[1];
rz(-1.5685969) q[3];
sx q[3];
rz(-1.723395) q[3];
sx q[3];
rz(-0.27657498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6039383) q[2];
sx q[2];
rz(-2.5014169) q[2];
sx q[2];
rz(-1.630265) q[2];
rz(2.9553735) q[3];
sx q[3];
rz(-2.7112466) q[3];
sx q[3];
rz(2.490624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6179825) q[0];
sx q[0];
rz(-1.1815434) q[0];
sx q[0];
rz(-3.004916) q[0];
rz(0.094820529) q[1];
sx q[1];
rz(-0.46351981) q[1];
sx q[1];
rz(-2.3725841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13327293) q[0];
sx q[0];
rz(-1.1964487) q[0];
sx q[0];
rz(-1.4248614) q[0];
x q[1];
rz(1.9273754) q[2];
sx q[2];
rz(-2.2327023) q[2];
sx q[2];
rz(1.8476768) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5040759) q[1];
sx q[1];
rz(-2.3055162) q[1];
sx q[1];
rz(-1.7903663) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3202928) q[3];
sx q[3];
rz(-3.0022096) q[3];
sx q[3];
rz(1.936862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7424221) q[2];
sx q[2];
rz(-2.6593282) q[2];
sx q[2];
rz(1.4313618) q[2];
rz(-2.6509905) q[3];
sx q[3];
rz(-2.6503745) q[3];
sx q[3];
rz(0.77559364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2630149) q[0];
sx q[0];
rz(-0.37913015) q[0];
sx q[0];
rz(2.0468792) q[0];
rz(-1.8393983) q[1];
sx q[1];
rz(-2.1523988) q[1];
sx q[1];
rz(-0.0040231752) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1510295) q[0];
sx q[0];
rz(-2.4044665) q[0];
sx q[0];
rz(-1.766979) q[0];
x q[1];
rz(-0.73765678) q[2];
sx q[2];
rz(-2.1263543) q[2];
sx q[2];
rz(0.45762024) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.04847542) q[1];
sx q[1];
rz(-2.8826536) q[1];
sx q[1];
rz(-0.11872752) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8588649) q[3];
sx q[3];
rz(-2.4077031) q[3];
sx q[3];
rz(-0.80724387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5502988) q[2];
sx q[2];
rz(-0.58612263) q[2];
sx q[2];
rz(2.6934521) q[2];
rz(2.6552933) q[3];
sx q[3];
rz(-2.2794162) q[3];
sx q[3];
rz(-1.1264616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22911856) q[0];
sx q[0];
rz(-1.4294701) q[0];
sx q[0];
rz(-2.4718156) q[0];
rz(1.229333) q[1];
sx q[1];
rz(-2.6826617) q[1];
sx q[1];
rz(-2.5965447) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16825039) q[0];
sx q[0];
rz(-2.0297883) q[0];
sx q[0];
rz(-0.58990546) q[0];
rz(-pi) q[1];
rz(-0.9093693) q[2];
sx q[2];
rz(-1.0535568) q[2];
sx q[2];
rz(1.6632207) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5601756) q[1];
sx q[1];
rz(-0.22318527) q[1];
sx q[1];
rz(-1.2944004) q[1];
x q[2];
rz(0.46579428) q[3];
sx q[3];
rz(-2.287103) q[3];
sx q[3];
rz(2.4972625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.58686078) q[2];
sx q[2];
rz(-1.6394337) q[2];
sx q[2];
rz(0.16793212) q[2];
rz(1.933291) q[3];
sx q[3];
rz(-2.8390563) q[3];
sx q[3];
rz(-1.4471853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55030441) q[0];
sx q[0];
rz(-1.8809603) q[0];
sx q[0];
rz(0.0038797832) q[0];
rz(0.30715352) q[1];
sx q[1];
rz(-0.75633621) q[1];
sx q[1];
rz(-0.53877962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85941797) q[0];
sx q[0];
rz(-0.27056405) q[0];
sx q[0];
rz(-2.5859588) q[0];
x q[1];
rz(2.8454119) q[2];
sx q[2];
rz(-2.6950172) q[2];
sx q[2];
rz(-2.2889529) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5153577) q[1];
sx q[1];
rz(-1.9222676) q[1];
sx q[1];
rz(2.5280747) q[1];
rz(-pi) q[2];
rz(2.410048) q[3];
sx q[3];
rz(-1.9172102) q[3];
sx q[3];
rz(-2.9820095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4469845) q[2];
sx q[2];
rz(-1.0993404) q[2];
sx q[2];
rz(1.2896607) q[2];
rz(-1.6192216) q[3];
sx q[3];
rz(-2.3963942) q[3];
sx q[3];
rz(1.4815909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8412142) q[0];
sx q[0];
rz(-0.9391681) q[0];
sx q[0];
rz(-3.0389431) q[0];
rz(0.73219055) q[1];
sx q[1];
rz(-0.79473549) q[1];
sx q[1];
rz(0.57714677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6711376) q[0];
sx q[0];
rz(-1.5744493) q[0];
sx q[0];
rz(-1.5994344) q[0];
rz(-pi) q[1];
rz(2.8961298) q[2];
sx q[2];
rz(-2.5405014) q[2];
sx q[2];
rz(-0.18537755) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6070093) q[1];
sx q[1];
rz(-1.680169) q[1];
sx q[1];
rz(-1.4916199) q[1];
rz(-pi) q[2];
rz(-2.7957245) q[3];
sx q[3];
rz(-0.78639275) q[3];
sx q[3];
rz(-2.6707471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.131549) q[2];
sx q[2];
rz(-0.71655822) q[2];
sx q[2];
rz(-1.6183759) q[2];
rz(0.067954436) q[3];
sx q[3];
rz(-2.8165635) q[3];
sx q[3];
rz(0.84097356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0129358) q[0];
sx q[0];
rz(-1.034863) q[0];
sx q[0];
rz(-2.3701684) q[0];
rz(0.5386638) q[1];
sx q[1];
rz(-1.3464876) q[1];
sx q[1];
rz(1.0661941) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5504042) q[0];
sx q[0];
rz(-0.26445358) q[0];
sx q[0];
rz(-1.868684) q[0];
x q[1];
rz(-1.9692405) q[2];
sx q[2];
rz(-0.5486998) q[2];
sx q[2];
rz(1.0469701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9061588) q[1];
sx q[1];
rz(-1.7879491) q[1];
sx q[1];
rz(1.557334) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3049502) q[3];
sx q[3];
rz(-1.2472594) q[3];
sx q[3];
rz(-0.27184799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.44925877) q[2];
sx q[2];
rz(-0.80654311) q[2];
sx q[2];
rz(-2.6991357) q[2];
rz(-0.19065204) q[3];
sx q[3];
rz(-2.4176044) q[3];
sx q[3];
rz(-0.75907069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8964748) q[0];
sx q[0];
rz(-0.18558311) q[0];
sx q[0];
rz(1.8316487) q[0];
rz(2.7438121) q[1];
sx q[1];
rz(-0.82292992) q[1];
sx q[1];
rz(1.7479053) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0017477) q[0];
sx q[0];
rz(-2.4187517) q[0];
sx q[0];
rz(-0.87840898) q[0];
x q[1];
rz(0.23598052) q[2];
sx q[2];
rz(-1.2771291) q[2];
sx q[2];
rz(0.27494173) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.85227784) q[1];
sx q[1];
rz(-2.5767269) q[1];
sx q[1];
rz(-0.58128618) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0354064) q[3];
sx q[3];
rz(-1.3732037) q[3];
sx q[3];
rz(1.1923517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.028367793) q[2];
sx q[2];
rz(-0.71030474) q[2];
sx q[2];
rz(-3.0681211) q[2];
rz(-0.69463378) q[3];
sx q[3];
rz(-1.8690542) q[3];
sx q[3];
rz(2.1139483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6756814) q[0];
sx q[0];
rz(-2.9407192) q[0];
sx q[0];
rz(2.1821816) q[0];
rz(2.361946) q[1];
sx q[1];
rz(-2.1387073) q[1];
sx q[1];
rz(-0.27228212) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1355541) q[0];
sx q[0];
rz(-1.557543) q[0];
sx q[0];
rz(-3.1363961) q[0];
rz(-pi) q[1];
rz(1.2113038) q[2];
sx q[2];
rz(-2.4036416) q[2];
sx q[2];
rz(2.8518845) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12676621) q[1];
sx q[1];
rz(-2.0527168) q[1];
sx q[1];
rz(-2.6792106) q[1];
rz(0.63377728) q[3];
sx q[3];
rz(-2.2867775) q[3];
sx q[3];
rz(-2.8785126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4043364) q[2];
sx q[2];
rz(-1.6360838) q[2];
sx q[2];
rz(-2.9340202) q[2];
rz(-0.10457822) q[3];
sx q[3];
rz(-0.50555491) q[3];
sx q[3];
rz(-1.526621) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4614457) q[0];
sx q[0];
rz(-1.2274281) q[0];
sx q[0];
rz(-2.0848059) q[0];
rz(0.54569221) q[1];
sx q[1];
rz(-0.97815198) q[1];
sx q[1];
rz(0.32870764) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3378819) q[0];
sx q[0];
rz(-2.0265371) q[0];
sx q[0];
rz(0.36949879) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4942371) q[2];
sx q[2];
rz(-3.050905) q[2];
sx q[2];
rz(0.91435223) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2519319) q[1];
sx q[1];
rz(-1.9128107) q[1];
sx q[1];
rz(-3.1238848) q[1];
x q[2];
rz(-0.77158934) q[3];
sx q[3];
rz(-1.3719012) q[3];
sx q[3];
rz(2.1344913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9823965) q[2];
sx q[2];
rz(-3.0089162) q[2];
sx q[2];
rz(-1.1245493) q[2];
rz(-3.0647965) q[3];
sx q[3];
rz(-2.7087637) q[3];
sx q[3];
rz(-0.92397773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15976739) q[0];
sx q[0];
rz(-1.2564909) q[0];
sx q[0];
rz(1.5782574) q[0];
rz(0.81623296) q[1];
sx q[1];
rz(-1.4030133) q[1];
sx q[1];
rz(-1.0907008) q[1];
rz(-2.0682206) q[2];
sx q[2];
rz(-2.2523027) q[2];
sx q[2];
rz(-0.61423617) q[2];
rz(-2.7362551) q[3];
sx q[3];
rz(-1.5673076) q[3];
sx q[3];
rz(0.81893541) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
