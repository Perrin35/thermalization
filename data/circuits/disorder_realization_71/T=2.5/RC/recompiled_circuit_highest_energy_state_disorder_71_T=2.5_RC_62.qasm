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
rz(0.55819297) q[0];
sx q[0];
rz(3.8605122) q[0];
sx q[0];
rz(10.100848) q[0];
rz(-2.8683635) q[1];
sx q[1];
rz(-0.9571119) q[1];
sx q[1];
rz(-0.91125429) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.229743) q[0];
sx q[0];
rz(-1.5485244) q[0];
sx q[0];
rz(-0.020072083) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5780294) q[2];
sx q[2];
rz(-0.94524318) q[2];
sx q[2];
rz(0.48239732) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.82750597) q[1];
sx q[1];
rz(-2.1438399) q[1];
sx q[1];
rz(-1.5324462) q[1];
rz(3.0082672) q[3];
sx q[3];
rz(-0.76078868) q[3];
sx q[3];
rz(-0.68851346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.91997826) q[2];
sx q[2];
rz(-2.6875434) q[2];
sx q[2];
rz(0.81217074) q[2];
rz(0.075856097) q[3];
sx q[3];
rz(-2.4935738) q[3];
sx q[3];
rz(0.59504741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6512063) q[0];
sx q[0];
rz(-2.6631329) q[0];
sx q[0];
rz(-1.2745717) q[0];
rz(-1.6365341) q[1];
sx q[1];
rz(-0.36205629) q[1];
sx q[1];
rz(1.1638181) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464435) q[0];
sx q[0];
rz(-2.4730485) q[0];
sx q[0];
rz(-2.3666381) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1357434) q[2];
sx q[2];
rz(-1.1052624) q[2];
sx q[2];
rz(-2.7201544) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5608985) q[1];
sx q[1];
rz(-2.1836062) q[1];
sx q[1];
rz(1.3152907) q[1];
rz(2.1970956) q[3];
sx q[3];
rz(-0.5118733) q[3];
sx q[3];
rz(-2.8688489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85094467) q[2];
sx q[2];
rz(-3.1182351) q[2];
sx q[2];
rz(-1.5604431) q[2];
rz(-0.23049878) q[3];
sx q[3];
rz(-1.0483619) q[3];
sx q[3];
rz(0.44052625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59349638) q[0];
sx q[0];
rz(-0.31262147) q[0];
sx q[0];
rz(0.12259677) q[0];
rz(2.3444046) q[1];
sx q[1];
rz(-2.1295348) q[1];
sx q[1];
rz(-0.0880934) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6428292) q[0];
sx q[0];
rz(-0.82422148) q[0];
sx q[0];
rz(2.7902664) q[0];
rz(3.0584774) q[2];
sx q[2];
rz(-1.3440455) q[2];
sx q[2];
rz(2.6672358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3077724) q[1];
sx q[1];
rz(-1.4188179) q[1];
sx q[1];
rz(0.13843468) q[1];
x q[2];
rz(-2.5949305) q[3];
sx q[3];
rz(-0.61184363) q[3];
sx q[3];
rz(0.33239588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.88283551) q[2];
sx q[2];
rz(-1.5828524) q[2];
sx q[2];
rz(1.3239512) q[2];
rz(-2.6435408) q[3];
sx q[3];
rz(-2.1784454) q[3];
sx q[3];
rz(0.74670416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8525456) q[0];
sx q[0];
rz(-0.15327029) q[0];
sx q[0];
rz(-2.9943941) q[0];
rz(2.4920801) q[1];
sx q[1];
rz(-1.5107061) q[1];
sx q[1];
rz(1.46896) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8627074) q[0];
sx q[0];
rz(-1.577958) q[0];
sx q[0];
rz(1.347752) q[0];
x q[1];
rz(-2.9584453) q[2];
sx q[2];
rz(-1.9948655) q[2];
sx q[2];
rz(1.1893502) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8007954) q[1];
sx q[1];
rz(-0.95424517) q[1];
sx q[1];
rz(-0.430213) q[1];
rz(-pi) q[2];
rz(0.62375237) q[3];
sx q[3];
rz(-2.459708) q[3];
sx q[3];
rz(0.15450134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2088251) q[2];
sx q[2];
rz(-2.6471477) q[2];
sx q[2];
rz(-0.24173582) q[2];
rz(2.406481) q[3];
sx q[3];
rz(-1.7416411) q[3];
sx q[3];
rz(-0.66489768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6786137) q[0];
sx q[0];
rz(-2.0942056) q[0];
sx q[0];
rz(-3.0188634) q[0];
rz(0.8853451) q[1];
sx q[1];
rz(-2.131999) q[1];
sx q[1];
rz(-2.5811894) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80622506) q[0];
sx q[0];
rz(-2.0809552) q[0];
sx q[0];
rz(2.9377743) q[0];
rz(-0.63626115) q[2];
sx q[2];
rz(-0.36379746) q[2];
sx q[2];
rz(-0.76278245) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97241114) q[1];
sx q[1];
rz(-1.3126846) q[1];
sx q[1];
rz(2.045579) q[1];
x q[2];
rz(2.6953648) q[3];
sx q[3];
rz(-1.130502) q[3];
sx q[3];
rz(-1.4780731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5237328) q[2];
sx q[2];
rz(-0.79404074) q[2];
sx q[2];
rz(-0.20475556) q[2];
rz(-2.2023885) q[3];
sx q[3];
rz(-1.3699646) q[3];
sx q[3];
rz(-0.088976629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.8646249) q[0];
sx q[0];
rz(-3.0245916) q[0];
sx q[0];
rz(0.65698874) q[0];
rz(2.9981546) q[1];
sx q[1];
rz(-1.5851494) q[1];
sx q[1];
rz(-0.73854804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8764832) q[0];
sx q[0];
rz(-0.68317693) q[0];
sx q[0];
rz(2.9580412) q[0];
rz(-0.22566585) q[2];
sx q[2];
rz(-1.0899001) q[2];
sx q[2];
rz(-1.4165598) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.025102928) q[1];
sx q[1];
rz(-1.2560351) q[1];
sx q[1];
rz(0.88309755) q[1];
rz(0.62885999) q[3];
sx q[3];
rz(-0.81506461) q[3];
sx q[3];
rz(3.0185543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2584381) q[2];
sx q[2];
rz(-2.482735) q[2];
sx q[2];
rz(-3.1385885) q[2];
rz(-2.352412) q[3];
sx q[3];
rz(-2.7572258) q[3];
sx q[3];
rz(1.9743617) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0810735) q[0];
sx q[0];
rz(-0.34927148) q[0];
sx q[0];
rz(-2.9935167) q[0];
rz(-0.5332467) q[1];
sx q[1];
rz(-0.24594578) q[1];
sx q[1];
rz(-1.6291133) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027372885) q[0];
sx q[0];
rz(-2.5974869) q[0];
sx q[0];
rz(-2.6065488) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.523411) q[2];
sx q[2];
rz(-0.91518171) q[2];
sx q[2];
rz(2.5822292) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.032669205) q[1];
sx q[1];
rz(-1.2315688) q[1];
sx q[1];
rz(0.12813385) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84866555) q[3];
sx q[3];
rz(-1.7300528) q[3];
sx q[3];
rz(-1.4383663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.44275722) q[2];
sx q[2];
rz(-1.4407225) q[2];
sx q[2];
rz(-3.0485145) q[2];
rz(-3.0794411) q[3];
sx q[3];
rz(-2.744894) q[3];
sx q[3];
rz(1.9183581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1208948) q[0];
sx q[0];
rz(-0.59497213) q[0];
sx q[0];
rz(2.4769532) q[0];
rz(1.7918034) q[1];
sx q[1];
rz(-2.2638958) q[1];
sx q[1];
rz(0.68914366) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6038325) q[0];
sx q[0];
rz(-1.6869154) q[0];
sx q[0];
rz(0.24531792) q[0];
rz(-2.2623863) q[2];
sx q[2];
rz(-1.9636646) q[2];
sx q[2];
rz(-0.61516064) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2011678) q[1];
sx q[1];
rz(-1.42527) q[1];
sx q[1];
rz(0.31647233) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5102622) q[3];
sx q[3];
rz(-1.406296) q[3];
sx q[3];
rz(-1.1806844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.050921116) q[2];
sx q[2];
rz(-0.65616578) q[2];
sx q[2];
rz(-1.3760759) q[2];
rz(-1.225166) q[3];
sx q[3];
rz(-2.4296032) q[3];
sx q[3];
rz(3.0998949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10251481) q[0];
sx q[0];
rz(-0.044487655) q[0];
sx q[0];
rz(0.59120375) q[0];
rz(-0.2419596) q[1];
sx q[1];
rz(-2.0457025) q[1];
sx q[1];
rz(2.718149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71899881) q[0];
sx q[0];
rz(-2.7072577) q[0];
sx q[0];
rz(2.3584764) q[0];
rz(-pi) q[1];
rz(0.4884004) q[2];
sx q[2];
rz(-1.6624221) q[2];
sx q[2];
rz(0.56885251) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8264173) q[1];
sx q[1];
rz(-1.2741544) q[1];
sx q[1];
rz(-2.1267932) q[1];
rz(2.5980392) q[3];
sx q[3];
rz(-0.69216903) q[3];
sx q[3];
rz(-0.33351937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5847136) q[2];
sx q[2];
rz(-2.5355279) q[2];
sx q[2];
rz(1.9752183) q[2];
rz(1.1276468) q[3];
sx q[3];
rz(-2.1731264) q[3];
sx q[3];
rz(-2.6166272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.128085) q[0];
sx q[0];
rz(-2.7078244) q[0];
sx q[0];
rz(0.67310131) q[0];
rz(2.290944) q[1];
sx q[1];
rz(-1.7885845) q[1];
sx q[1];
rz(-2.8180715) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0304383) q[0];
sx q[0];
rz(-1.1980828) q[0];
sx q[0];
rz(-0.26813676) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8658429) q[2];
sx q[2];
rz(-2.033503) q[2];
sx q[2];
rz(1.5327061) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4648351) q[1];
sx q[1];
rz(-0.82993648) q[1];
sx q[1];
rz(-2.6840032) q[1];
x q[2];
rz(3.1296842) q[3];
sx q[3];
rz(-0.58387305) q[3];
sx q[3];
rz(-1.5990822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9234151) q[2];
sx q[2];
rz(-2.9604993) q[2];
sx q[2];
rz(0.072048135) q[2];
rz(2.1273023) q[3];
sx q[3];
rz(-0.91599661) q[3];
sx q[3];
rz(-2.5924276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88076787) q[0];
sx q[0];
rz(-1.4243955) q[0];
sx q[0];
rz(-2.0212174) q[0];
rz(2.8063759) q[1];
sx q[1];
rz(-2.4991279) q[1];
sx q[1];
rz(2.0229708) q[1];
rz(-1.3041244) q[2];
sx q[2];
rz(-0.49378569) q[2];
sx q[2];
rz(2.7616382) q[2];
rz(1.5961965) q[3];
sx q[3];
rz(-1.7138154) q[3];
sx q[3];
rz(-2.9830767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
