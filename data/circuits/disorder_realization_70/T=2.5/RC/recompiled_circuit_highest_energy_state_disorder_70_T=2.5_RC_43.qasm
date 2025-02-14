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
rz(-0.47082585) q[0];
sx q[0];
rz(-2.6915221) q[0];
sx q[0];
rz(0.39029628) q[0];
rz(0.14416873) q[1];
sx q[1];
rz(-1.498797) q[1];
sx q[1];
rz(2.1013451) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1821572) q[0];
sx q[0];
rz(-2.6819026) q[0];
sx q[0];
rz(0.88031405) q[0];
x q[1];
rz(0.36379142) q[2];
sx q[2];
rz(-1.7161088) q[2];
sx q[2];
rz(2.7940968) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92683243) q[1];
sx q[1];
rz(-0.70509395) q[1];
sx q[1];
rz(1.2486723) q[1];
rz(-1.5055429) q[3];
sx q[3];
rz(-2.3480519) q[3];
sx q[3];
rz(-0.82543594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.9425923) q[2];
sx q[2];
rz(-2.5358443) q[2];
sx q[2];
rz(1.595363) q[2];
rz(-0.47510251) q[3];
sx q[3];
rz(-0.64517704) q[3];
sx q[3];
rz(1.9861541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5763181) q[0];
sx q[0];
rz(-2.1529038) q[0];
sx q[0];
rz(-0.43462547) q[0];
rz(-1.0644396) q[1];
sx q[1];
rz(-0.39355215) q[1];
sx q[1];
rz(-0.89148608) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8572509) q[0];
sx q[0];
rz(-1.9280757) q[0];
sx q[0];
rz(2.7140359) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72618809) q[2];
sx q[2];
rz(-1.9627769) q[2];
sx q[2];
rz(-1.7418944) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3825917) q[1];
sx q[1];
rz(-1.0487952) q[1];
sx q[1];
rz(-2.0981789) q[1];
rz(-pi) q[2];
rz(0.029835506) q[3];
sx q[3];
rz(-2.2298457) q[3];
sx q[3];
rz(-1.713879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.50856227) q[2];
sx q[2];
rz(-0.56266251) q[2];
sx q[2];
rz(1.0594581) q[2];
rz(1.2262454) q[3];
sx q[3];
rz(-1.6112593) q[3];
sx q[3];
rz(-3.0219769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0480963) q[0];
sx q[0];
rz(-2.7432848) q[0];
sx q[0];
rz(-0.68921971) q[0];
rz(-2.7293909) q[1];
sx q[1];
rz(-2.0239794) q[1];
sx q[1];
rz(-1.9788007) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49814082) q[0];
sx q[0];
rz(-1.2664794) q[0];
sx q[0];
rz(2.9813719) q[0];
rz(-pi) q[1];
rz(-2.3648889) q[2];
sx q[2];
rz(-1.113184) q[2];
sx q[2];
rz(2.6991778) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38123576) q[1];
sx q[1];
rz(-2.3489526) q[1];
sx q[1];
rz(2.4355678) q[1];
rz(-pi) q[2];
rz(2.8619295) q[3];
sx q[3];
rz(-1.8496017) q[3];
sx q[3];
rz(-1.6300843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73843655) q[2];
sx q[2];
rz(-1.7313892) q[2];
sx q[2];
rz(-0.23846826) q[2];
rz(2.9078935) q[3];
sx q[3];
rz(-0.77478474) q[3];
sx q[3];
rz(0.15271798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9018263) q[0];
sx q[0];
rz(-0.16401839) q[0];
sx q[0];
rz(-2.8357586) q[0];
rz(0.26890525) q[1];
sx q[1];
rz(-2.3482359) q[1];
sx q[1];
rz(-2.7893524) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7765609) q[0];
sx q[0];
rz(-1.7045341) q[0];
sx q[0];
rz(3.1015009) q[0];
rz(-pi) q[1];
rz(2.1969123) q[2];
sx q[2];
rz(-1.4123823) q[2];
sx q[2];
rz(0.84412727) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8246987) q[1];
sx q[1];
rz(-2.3815063) q[1];
sx q[1];
rz(-0.58652564) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8958578) q[3];
sx q[3];
rz(-0.80348368) q[3];
sx q[3];
rz(-1.5218228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8504146) q[2];
sx q[2];
rz(-1.520227) q[2];
sx q[2];
rz(-2.9746383) q[2];
rz(3.0252365) q[3];
sx q[3];
rz(-0.51893026) q[3];
sx q[3];
rz(-1.8714582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6312234) q[0];
sx q[0];
rz(-0.33741697) q[0];
sx q[0];
rz(1.1153197) q[0];
rz(1.2315617) q[1];
sx q[1];
rz(-1.7111338) q[1];
sx q[1];
rz(-0.11428741) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63681215) q[0];
sx q[0];
rz(-2.0858602) q[0];
sx q[0];
rz(-2.5367141) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27844001) q[2];
sx q[2];
rz(-1.2781218) q[2];
sx q[2];
rz(-0.038624374) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48620079) q[1];
sx q[1];
rz(-0.95410141) q[1];
sx q[1];
rz(-2.0808897) q[1];
rz(-pi) q[2];
rz(-1.565236) q[3];
sx q[3];
rz(-0.86478028) q[3];
sx q[3];
rz(1.8354285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.085122434) q[2];
sx q[2];
rz(-1.9279927) q[2];
sx q[2];
rz(-1.6950133) q[2];
rz(2.1794686) q[3];
sx q[3];
rz(-0.4947997) q[3];
sx q[3];
rz(-1.3303293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98244786) q[0];
sx q[0];
rz(-1.5788989) q[0];
sx q[0];
rz(-2.6084117) q[0];
rz(1.606733) q[1];
sx q[1];
rz(-2.3695562) q[1];
sx q[1];
rz(0.74347043) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0395775) q[0];
sx q[0];
rz(-1.8626889) q[0];
sx q[0];
rz(1.6534228) q[0];
rz(1.1986046) q[2];
sx q[2];
rz(-0.57665885) q[2];
sx q[2];
rz(-1.8221591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4308692) q[1];
sx q[1];
rz(-0.9473045) q[1];
sx q[1];
rz(1.2602978) q[1];
x q[2];
rz(-0.95249259) q[3];
sx q[3];
rz(-1.8466766) q[3];
sx q[3];
rz(-0.80989686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6712436) q[2];
sx q[2];
rz(-0.84731421) q[2];
sx q[2];
rz(-2.642855) q[2];
rz(0.70164743) q[3];
sx q[3];
rz(-2.4528153) q[3];
sx q[3];
rz(-0.65830314) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43211234) q[0];
sx q[0];
rz(-1.226959) q[0];
sx q[0];
rz(-3.0416601) q[0];
rz(-1.8607032) q[1];
sx q[1];
rz(-0.51191267) q[1];
sx q[1];
rz(-2.4023712) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4866132) q[0];
sx q[0];
rz(-0.57679048) q[0];
sx q[0];
rz(-0.77107112) q[0];
rz(-0.03050892) q[2];
sx q[2];
rz(-2.0000474) q[2];
sx q[2];
rz(-0.78684083) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2120253) q[1];
sx q[1];
rz(-1.9627357) q[1];
sx q[1];
rz(2.2076616) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5600863) q[3];
sx q[3];
rz(-2.3751343) q[3];
sx q[3];
rz(-0.86957726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7382536) q[2];
sx q[2];
rz(-0.6311987) q[2];
sx q[2];
rz(-1.6048019) q[2];
rz(-1.4295476) q[3];
sx q[3];
rz(-1.6481383) q[3];
sx q[3];
rz(2.3454989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1256063) q[0];
sx q[0];
rz(-0.37785372) q[0];
sx q[0];
rz(2.0269537) q[0];
rz(2.5573225) q[1];
sx q[1];
rz(-0.92002267) q[1];
sx q[1];
rz(3.027473) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18543359) q[0];
sx q[0];
rz(-1.2853649) q[0];
sx q[0];
rz(0.6298084) q[0];
rz(-pi) q[1];
rz(-1.0074963) q[2];
sx q[2];
rz(-0.85342583) q[2];
sx q[2];
rz(-0.12441758) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89592592) q[1];
sx q[1];
rz(-1.1628431) q[1];
sx q[1];
rz(-2.308564) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0242943) q[3];
sx q[3];
rz(-1.3051194) q[3];
sx q[3];
rz(1.0289711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.70322651) q[2];
sx q[2];
rz(-0.72202903) q[2];
sx q[2];
rz(-0.88877338) q[2];
rz(-1.5416386) q[3];
sx q[3];
rz(-1.2069353) q[3];
sx q[3];
rz(0.82109872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.090076598) q[0];
sx q[0];
rz(-0.30802825) q[0];
sx q[0];
rz(-2.8544881) q[0];
rz(1.9376532) q[1];
sx q[1];
rz(-2.2724889) q[1];
sx q[1];
rz(-1.513011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2503165) q[0];
sx q[0];
rz(-2.7990632) q[0];
sx q[0];
rz(-2.9129759) q[0];
x q[1];
rz(-2.2193935) q[2];
sx q[2];
rz(-1.2051799) q[2];
sx q[2];
rz(0.43482414) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2681458) q[1];
sx q[1];
rz(-1.526274) q[1];
sx q[1];
rz(-2.4067307) q[1];
rz(-pi) q[2];
rz(-0.00079454409) q[3];
sx q[3];
rz(-1.5738788) q[3];
sx q[3];
rz(2.2828988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63605753) q[2];
sx q[2];
rz(-1.1670185) q[2];
sx q[2];
rz(-2.8324845) q[2];
rz(1.9348034) q[3];
sx q[3];
rz(-0.38201067) q[3];
sx q[3];
rz(2.8822854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9647656) q[0];
sx q[0];
rz(-0.72787705) q[0];
sx q[0];
rz(2.4859909) q[0];
rz(0.62878311) q[1];
sx q[1];
rz(-1.4097593) q[1];
sx q[1];
rz(-1.383673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1882598) q[0];
sx q[0];
rz(-1.6142705) q[0];
sx q[0];
rz(-0.047634634) q[0];
rz(-pi) q[1];
rz(-0.45510095) q[2];
sx q[2];
rz(-0.88717194) q[2];
sx q[2];
rz(3.07651) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1392581) q[1];
sx q[1];
rz(-2.9911925) q[1];
sx q[1];
rz(-0.12087442) q[1];
rz(-pi) q[2];
rz(-2.6419956) q[3];
sx q[3];
rz(-1.4471421) q[3];
sx q[3];
rz(-0.62859231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3722966) q[2];
sx q[2];
rz(-1.5089704) q[2];
sx q[2];
rz(0.24202913) q[2];
rz(0.58487839) q[3];
sx q[3];
rz(-2.4681028) q[3];
sx q[3];
rz(2.6175595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.723421) q[0];
sx q[0];
rz(-0.62234288) q[0];
sx q[0];
rz(-0.55707669) q[0];
rz(0.88049018) q[1];
sx q[1];
rz(-1.7387895) q[1];
sx q[1];
rz(-2.1008076) q[1];
rz(1.6014848) q[2];
sx q[2];
rz(-2.9353113) q[2];
sx q[2];
rz(-2.7842709) q[2];
rz(-1.7193033) q[3];
sx q[3];
rz(-2.5703493) q[3];
sx q[3];
rz(-0.66019365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
