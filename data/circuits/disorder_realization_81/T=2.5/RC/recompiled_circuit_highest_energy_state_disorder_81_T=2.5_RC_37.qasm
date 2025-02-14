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
rz(2.7414275) q[0];
sx q[0];
rz(-2.6097809) q[0];
sx q[0];
rz(-2.3398633) q[0];
rz(-2.5126558) q[1];
sx q[1];
rz(-1.3087654) q[1];
sx q[1];
rz(-0.3892678) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67943924) q[0];
sx q[0];
rz(-0.91374699) q[0];
sx q[0];
rz(1.2253967) q[0];
rz(-2.9603329) q[2];
sx q[2];
rz(-1.6158293) q[2];
sx q[2];
rz(-0.35895106) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3013346) q[1];
sx q[1];
rz(-1.0691027) q[1];
sx q[1];
rz(-2.4444716) q[1];
rz(-2.2981529) q[3];
sx q[3];
rz(-1.6004082) q[3];
sx q[3];
rz(-0.51472291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6486711) q[2];
sx q[2];
rz(-0.45009437) q[2];
sx q[2];
rz(0.2613124) q[2];
rz(2.785545) q[3];
sx q[3];
rz(-1.1084403) q[3];
sx q[3];
rz(2.0480919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68591958) q[0];
sx q[0];
rz(-2.5643667) q[0];
sx q[0];
rz(-2.8976231) q[0];
rz(2.7120554) q[1];
sx q[1];
rz(-1.669084) q[1];
sx q[1];
rz(-2.4006749) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2527494) q[0];
sx q[0];
rz(-2.6359632) q[0];
sx q[0];
rz(-2.3028785) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6119611) q[2];
sx q[2];
rz(-1.1897693) q[2];
sx q[2];
rz(1.6356026) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.59830211) q[1];
sx q[1];
rz(-1.0183698) q[1];
sx q[1];
rz(-1.8364947) q[1];
rz(-pi) q[2];
rz(0.55161019) q[3];
sx q[3];
rz(-0.79539652) q[3];
sx q[3];
rz(-1.7097527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4247158) q[2];
sx q[2];
rz(-1.6191142) q[2];
sx q[2];
rz(0.87872046) q[2];
rz(-0.97484136) q[3];
sx q[3];
rz(-1.2474371) q[3];
sx q[3];
rz(-1.6802906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9413302) q[0];
sx q[0];
rz(-2.5859079) q[0];
sx q[0];
rz(2.9479807) q[0];
rz(-3.0337453) q[1];
sx q[1];
rz(-0.71433181) q[1];
sx q[1];
rz(-2.329619) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40859544) q[0];
sx q[0];
rz(-0.46950668) q[0];
sx q[0];
rz(1.3256737) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0103364) q[2];
sx q[2];
rz(-2.3046846) q[2];
sx q[2];
rz(0.79272717) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4311422) q[1];
sx q[1];
rz(-2.6796088) q[1];
sx q[1];
rz(-1.0313018) q[1];
rz(-pi) q[2];
rz(-2.009684) q[3];
sx q[3];
rz(-1.697007) q[3];
sx q[3];
rz(0.47078339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0783483) q[2];
sx q[2];
rz(-1.8599267) q[2];
sx q[2];
rz(3.1380624) q[2];
rz(2.6053536) q[3];
sx q[3];
rz(-0.25501525) q[3];
sx q[3];
rz(-0.53853881) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2980767) q[0];
sx q[0];
rz(-1.902782) q[0];
sx q[0];
rz(-2.3115944) q[0];
rz(1.7371381) q[1];
sx q[1];
rz(-0.91157118) q[1];
sx q[1];
rz(-2.8703168) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0381841) q[0];
sx q[0];
rz(-1.8233946) q[0];
sx q[0];
rz(-0.20471666) q[0];
x q[1];
rz(0.6649632) q[2];
sx q[2];
rz(-1.5100484) q[2];
sx q[2];
rz(0.052968979) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2597166) q[1];
sx q[1];
rz(-1.9278942) q[1];
sx q[1];
rz(1.5413589) q[1];
x q[2];
rz(-2.3448571) q[3];
sx q[3];
rz(-1.6816144) q[3];
sx q[3];
rz(-3.022242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6259049) q[2];
sx q[2];
rz(-0.38334623) q[2];
sx q[2];
rz(2.3670727) q[2];
rz(-1.0033311) q[3];
sx q[3];
rz(-2.7480875) q[3];
sx q[3];
rz(0.5459319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14458732) q[0];
sx q[0];
rz(-2.044675) q[0];
sx q[0];
rz(2.0487832) q[0];
rz(-1.8057168) q[1];
sx q[1];
rz(-0.74936682) q[1];
sx q[1];
rz(-0.54378477) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0896801) q[0];
sx q[0];
rz(-2.8767715) q[0];
sx q[0];
rz(-1.7811243) q[0];
rz(0.67853482) q[2];
sx q[2];
rz(-1.6521875) q[2];
sx q[2];
rz(-2.396109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.005086686) q[1];
sx q[1];
rz(-0.6906116) q[1];
sx q[1];
rz(2.4588799) q[1];
x q[2];
rz(1.5690992) q[3];
sx q[3];
rz(-2.5212396) q[3];
sx q[3];
rz(1.2173295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.9076697) q[2];
sx q[2];
rz(-0.69207865) q[2];
sx q[2];
rz(2.862759) q[2];
rz(-0.28576609) q[3];
sx q[3];
rz(-2.2279492) q[3];
sx q[3];
rz(0.41551503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7061507) q[0];
sx q[0];
rz(-0.64318648) q[0];
sx q[0];
rz(1.0253133) q[0];
rz(3.0030491) q[1];
sx q[1];
rz(-1.5516611) q[1];
sx q[1];
rz(2.9655546) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6556102) q[0];
sx q[0];
rz(-1.8261251) q[0];
sx q[0];
rz(0.53389733) q[0];
rz(1.7654714) q[2];
sx q[2];
rz(-1.9241714) q[2];
sx q[2];
rz(-0.88549685) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61689395) q[1];
sx q[1];
rz(-2.0512676) q[1];
sx q[1];
rz(-2.8887649) q[1];
x q[2];
rz(-1.995924) q[3];
sx q[3];
rz(-1.7264724) q[3];
sx q[3];
rz(2.3666414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90870086) q[2];
sx q[2];
rz(-0.42753926) q[2];
sx q[2];
rz(-0.78988451) q[2];
rz(2.8120153) q[3];
sx q[3];
rz(-1.3720082) q[3];
sx q[3];
rz(-0.65143877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4306507) q[0];
sx q[0];
rz(-0.21345226) q[0];
sx q[0];
rz(2.1479204) q[0];
rz(1.0020533) q[1];
sx q[1];
rz(-0.99169815) q[1];
sx q[1];
rz(-1.6755982) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1260557) q[0];
sx q[0];
rz(-1.7483341) q[0];
sx q[0];
rz(2.4734495) q[0];
x q[1];
rz(-3.0266031) q[2];
sx q[2];
rz(-1.1144979) q[2];
sx q[2];
rz(1.4062238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42211275) q[1];
sx q[1];
rz(-1.7168448) q[1];
sx q[1];
rz(-2.9587347) q[1];
rz(-pi) q[2];
rz(-0.44963533) q[3];
sx q[3];
rz(-0.88346601) q[3];
sx q[3];
rz(0.80865142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7953636) q[2];
sx q[2];
rz(-0.84285223) q[2];
sx q[2];
rz(2.8818434) q[2];
rz(2.3527457) q[3];
sx q[3];
rz(-0.84281054) q[3];
sx q[3];
rz(-1.747725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2271093) q[0];
sx q[0];
rz(-1.8552584) q[0];
sx q[0];
rz(-0.9911384) q[0];
rz(-1.8046851) q[1];
sx q[1];
rz(-0.82695812) q[1];
sx q[1];
rz(1.7327259) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73118657) q[0];
sx q[0];
rz(-2.0265739) q[0];
sx q[0];
rz(-1.6996933) q[0];
rz(-pi) q[1];
rz(-2.1091923) q[2];
sx q[2];
rz(-1.4941895) q[2];
sx q[2];
rz(-1.3860764) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.046890586) q[1];
sx q[1];
rz(-0.45744236) q[1];
sx q[1];
rz(2.4777458) q[1];
rz(-pi) q[2];
rz(1.7400319) q[3];
sx q[3];
rz(-0.90247655) q[3];
sx q[3];
rz(3.075222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9764497) q[2];
sx q[2];
rz(-1.6565448) q[2];
sx q[2];
rz(0.58839798) q[2];
rz(0.31383651) q[3];
sx q[3];
rz(-2.2321759) q[3];
sx q[3];
rz(-2.5294144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2345851) q[0];
sx q[0];
rz(-1.122965) q[0];
sx q[0];
rz(-3.0871952) q[0];
rz(2.3811978) q[1];
sx q[1];
rz(-2.1957928) q[1];
sx q[1];
rz(-2.0423582) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71610036) q[0];
sx q[0];
rz(-1.4045225) q[0];
sx q[0];
rz(0.99207741) q[0];
x q[1];
rz(0.77950355) q[2];
sx q[2];
rz(-2.1981249) q[2];
sx q[2];
rz(1.9169538) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9304946) q[1];
sx q[1];
rz(-1.7451442) q[1];
sx q[1];
rz(-2.3424108) q[1];
rz(0.051088574) q[3];
sx q[3];
rz(-1.4600442) q[3];
sx q[3];
rz(0.95789078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86193639) q[2];
sx q[2];
rz(-1.6013689) q[2];
sx q[2];
rz(2.3242365) q[2];
rz(-0.014184626) q[3];
sx q[3];
rz(-2.7458906) q[3];
sx q[3];
rz(2.0524249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.31187439) q[0];
sx q[0];
rz(-2.0901966) q[0];
sx q[0];
rz(0.22148393) q[0];
rz(-0.43137506) q[1];
sx q[1];
rz(-2.2384877) q[1];
sx q[1];
rz(-2.121714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021602623) q[0];
sx q[0];
rz(-0.051030397) q[0];
sx q[0];
rz(2.7330815) q[0];
rz(-1.6758317) q[2];
sx q[2];
rz(-1.7713462) q[2];
sx q[2];
rz(0.54947621) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74032593) q[1];
sx q[1];
rz(-1.6049084) q[1];
sx q[1];
rz(1.9032065) q[1];
rz(-pi) q[2];
rz(1.7027036) q[3];
sx q[3];
rz(-1.7057837) q[3];
sx q[3];
rz(-1.5582946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43883103) q[2];
sx q[2];
rz(-1.2063382) q[2];
sx q[2];
rz(1.8863401) q[2];
rz(0.025764763) q[3];
sx q[3];
rz(-1.3296826) q[3];
sx q[3];
rz(-0.51938957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7230566) q[0];
sx q[0];
rz(-0.78934712) q[0];
sx q[0];
rz(-2.017372) q[0];
rz(2.8070246) q[1];
sx q[1];
rz(-2.6523013) q[1];
sx q[1];
rz(2.1025067) q[1];
rz(2.6417021) q[2];
sx q[2];
rz(-0.51639787) q[2];
sx q[2];
rz(-2.7698386) q[2];
rz(0.48503899) q[3];
sx q[3];
rz(-2.1342604) q[3];
sx q[3];
rz(0.66935183) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
