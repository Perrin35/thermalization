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
rz(-2.7201535) q[0];
sx q[0];
rz(-1.8869737) q[0];
sx q[0];
rz(-0.52809554) q[0];
rz(-0.34673196) q[1];
sx q[1];
rz(-1.81387) q[1];
sx q[1];
rz(-0.2833856) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58812773) q[0];
sx q[0];
rz(-1.4736543) q[0];
sx q[0];
rz(2.5573822) q[0];
x q[1];
rz(-2.0714226) q[2];
sx q[2];
rz(-2.3267022) q[2];
sx q[2];
rz(1.3484777) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8979075) q[1];
sx q[1];
rz(-2.7449611) q[1];
sx q[1];
rz(-0.70901386) q[1];
rz(-1.3259726) q[3];
sx q[3];
rz(-0.41327616) q[3];
sx q[3];
rz(-2.5101889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1957207) q[2];
sx q[2];
rz(-1.6678145) q[2];
sx q[2];
rz(-2.8742068) q[2];
rz(2.9684559) q[3];
sx q[3];
rz(-2.0045529) q[3];
sx q[3];
rz(-2.4587629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6807569) q[0];
sx q[0];
rz(-0.84592485) q[0];
sx q[0];
rz(2.8435006) q[0];
rz(0.2671034) q[1];
sx q[1];
rz(-0.88228455) q[1];
sx q[1];
rz(-0.90046901) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1698477) q[0];
sx q[0];
rz(-1.8119703) q[0];
sx q[0];
rz(1.3154047) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.995558) q[2];
sx q[2];
rz(-1.7196557) q[2];
sx q[2];
rz(1.7796041) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5067277) q[1];
sx q[1];
rz(-1.7657451) q[1];
sx q[1];
rz(-0.69201236) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9252842) q[3];
sx q[3];
rz(-2.1976314) q[3];
sx q[3];
rz(1.1594866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4233826) q[2];
sx q[2];
rz(-0.98576468) q[2];
sx q[2];
rz(-2.5441489) q[2];
rz(-1.0159703) q[3];
sx q[3];
rz(-1.6705325) q[3];
sx q[3];
rz(-1.5684675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628767) q[0];
sx q[0];
rz(-0.91575423) q[0];
sx q[0];
rz(0.065091982) q[0];
rz(2.0386631) q[1];
sx q[1];
rz(-1.256559) q[1];
sx q[1];
rz(2.7535313) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86923118) q[0];
sx q[0];
rz(-0.80565208) q[0];
sx q[0];
rz(0.3262354) q[0];
rz(1.4955639) q[2];
sx q[2];
rz(-1.2047775) q[2];
sx q[2];
rz(-0.96142069) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2414723) q[1];
sx q[1];
rz(-2.6992976) q[1];
sx q[1];
rz(-1.9458179) q[1];
rz(-pi) q[2];
rz(-0.89395071) q[3];
sx q[3];
rz(-0.36374796) q[3];
sx q[3];
rz(-1.8574024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1661561) q[2];
sx q[2];
rz(-2.1076951) q[2];
sx q[2];
rz(3.0530829) q[2];
rz(0.65781188) q[3];
sx q[3];
rz(-2.7031873) q[3];
sx q[3];
rz(-1.2737761) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45519644) q[0];
sx q[0];
rz(-0.96298591) q[0];
sx q[0];
rz(-0.92940593) q[0];
rz(1.9724253) q[1];
sx q[1];
rz(-2.2678352) q[1];
sx q[1];
rz(3.015231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0189782) q[0];
sx q[0];
rz(-0.41176957) q[0];
sx q[0];
rz(1.7420578) q[0];
x q[1];
rz(1.7772555) q[2];
sx q[2];
rz(-1.8935618) q[2];
sx q[2];
rz(-1.9545123) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.11495464) q[1];
sx q[1];
rz(-1.2628139) q[1];
sx q[1];
rz(2.6796209) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7533412) q[3];
sx q[3];
rz(-1.3352003) q[3];
sx q[3];
rz(0.65956957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.29752842) q[2];
sx q[2];
rz(-2.0487831) q[2];
sx q[2];
rz(-1.4275985) q[2];
rz(1.1045688) q[3];
sx q[3];
rz(-1.5449056) q[3];
sx q[3];
rz(0.67232084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4349058) q[0];
sx q[0];
rz(-2.4770985) q[0];
sx q[0];
rz(-1.8988761) q[0];
rz(-0.13936123) q[1];
sx q[1];
rz(-2.8262973) q[1];
sx q[1];
rz(-2.349966) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46027943) q[0];
sx q[0];
rz(-0.15817197) q[0];
sx q[0];
rz(1.3307443) q[0];
rz(-pi) q[1];
rz(-0.56812079) q[2];
sx q[2];
rz(-1.4213398) q[2];
sx q[2];
rz(0.40772256) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.086603786) q[1];
sx q[1];
rz(-1.5947002) q[1];
sx q[1];
rz(-1.2633282) q[1];
rz(0.11939998) q[3];
sx q[3];
rz(-0.97115937) q[3];
sx q[3];
rz(-1.7922082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41607729) q[2];
sx q[2];
rz(-1.6107586) q[2];
sx q[2];
rz(-1.3610972) q[2];
rz(2.1040037) q[3];
sx q[3];
rz(-1.6890084) q[3];
sx q[3];
rz(3.1328372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1182275) q[0];
sx q[0];
rz(-0.26015493) q[0];
sx q[0];
rz(0.17876974) q[0];
rz(0.44890064) q[1];
sx q[1];
rz(-1.9763549) q[1];
sx q[1];
rz(-1.9662205) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3942053) q[0];
sx q[0];
rz(-2.3763467) q[0];
sx q[0];
rz(3.0470362) q[0];
rz(1.825338) q[2];
sx q[2];
rz(-0.78331982) q[2];
sx q[2];
rz(-2.6412727) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0302802) q[1];
sx q[1];
rz(-0.86568173) q[1];
sx q[1];
rz(-2.6837803) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6244632) q[3];
sx q[3];
rz(-1.367336) q[3];
sx q[3];
rz(0.77282016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7281404) q[2];
sx q[2];
rz(-0.38273013) q[2];
sx q[2];
rz(-0.21391301) q[2];
rz(1.4296069) q[3];
sx q[3];
rz(-1.471328) q[3];
sx q[3];
rz(-1.3235929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0578617) q[0];
sx q[0];
rz(-1.2849176) q[0];
sx q[0];
rz(2.4833552) q[0];
rz(1.7103051) q[1];
sx q[1];
rz(-0.88600102) q[1];
sx q[1];
rz(-1.5586982) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3806136) q[0];
sx q[0];
rz(-1.7581994) q[0];
sx q[0];
rz(-0.44664573) q[0];
rz(-pi) q[1];
rz(0.6391292) q[2];
sx q[2];
rz(-0.37593146) q[2];
sx q[2];
rz(-1.3522918) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78608785) q[1];
sx q[1];
rz(-1.9078603) q[1];
sx q[1];
rz(-0.33871594) q[1];
rz(-pi) q[2];
rz(-2.57929) q[3];
sx q[3];
rz(-2.4829757) q[3];
sx q[3];
rz(-0.035619481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3805286) q[2];
sx q[2];
rz(-1.5917799) q[2];
sx q[2];
rz(-1.2065678) q[2];
rz(-1.0208463) q[3];
sx q[3];
rz(-2.6267093) q[3];
sx q[3];
rz(2.2339039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9970488) q[0];
sx q[0];
rz(-0.21268614) q[0];
sx q[0];
rz(-0.3120684) q[0];
rz(0.32876217) q[1];
sx q[1];
rz(-1.6203251) q[1];
sx q[1];
rz(-2.388248) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8327717) q[0];
sx q[0];
rz(-1.0139019) q[0];
sx q[0];
rz(2.4514319) q[0];
x q[1];
rz(-1.6095095) q[2];
sx q[2];
rz(-0.51100327) q[2];
sx q[2];
rz(-1.7450361) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4792013) q[1];
sx q[1];
rz(-2.0298784) q[1];
sx q[1];
rz(-2.4513112) q[1];
x q[2];
rz(-2.6238717) q[3];
sx q[3];
rz(-2.0690527) q[3];
sx q[3];
rz(2.7473828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3030777) q[2];
sx q[2];
rz(-1.7531771) q[2];
sx q[2];
rz(0.77084368) q[2];
rz(0.54909697) q[3];
sx q[3];
rz(-1.2884459) q[3];
sx q[3];
rz(-0.87636605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2751806) q[0];
sx q[0];
rz(-2.488945) q[0];
sx q[0];
rz(-1.3394248) q[0];
rz(-2.0088947) q[1];
sx q[1];
rz(-0.39342543) q[1];
sx q[1];
rz(3.0963617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68494481) q[0];
sx q[0];
rz(-0.47950011) q[0];
sx q[0];
rz(-1.142169) q[0];
rz(-pi) q[1];
x q[1];
rz(1.253183) q[2];
sx q[2];
rz(-0.3738974) q[2];
sx q[2];
rz(1.31391) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.69581) q[1];
sx q[1];
rz(-0.28430609) q[1];
sx q[1];
rz(2.6587756) q[1];
rz(2.8969403) q[3];
sx q[3];
rz(-2.1095697) q[3];
sx q[3];
rz(-0.95726171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2356147) q[2];
sx q[2];
rz(-2.2048042) q[2];
sx q[2];
rz(2.3255685) q[2];
rz(1.3927381) q[3];
sx q[3];
rz(-2.0413028) q[3];
sx q[3];
rz(1.6284774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29869646) q[0];
sx q[0];
rz(-0.5492292) q[0];
sx q[0];
rz(1.5811051) q[0];
rz(-0.16079482) q[1];
sx q[1];
rz(-2.000688) q[1];
sx q[1];
rz(-2.2198832) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.164249) q[0];
sx q[0];
rz(-1.5802339) q[0];
sx q[0];
rz(-1.759985) q[0];
rz(-pi) q[1];
rz(-2.5530898) q[2];
sx q[2];
rz(-1.0784618) q[2];
sx q[2];
rz(2.0534317) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8556577) q[1];
sx q[1];
rz(-1.8536356) q[1];
sx q[1];
rz(-2.9873579) q[1];
x q[2];
rz(0.21375938) q[3];
sx q[3];
rz(-1.4204657) q[3];
sx q[3];
rz(-1.1663471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2738652) q[2];
sx q[2];
rz(-2.2160857) q[2];
sx q[2];
rz(-2.0694536) q[2];
rz(-2.9169361) q[3];
sx q[3];
rz(-1.4916689) q[3];
sx q[3];
rz(0.91607654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0199725) q[0];
sx q[0];
rz(-2.2157123) q[0];
sx q[0];
rz(-2.8257688) q[0];
rz(0.074180457) q[1];
sx q[1];
rz(-1.0646432) q[1];
sx q[1];
rz(-0.021312996) q[1];
rz(-1.5431719) q[2];
sx q[2];
rz(-0.78747126) q[2];
sx q[2];
rz(-3.1046545) q[2];
rz(-2.3754121) q[3];
sx q[3];
rz(-1.2824367) q[3];
sx q[3];
rz(-2.5855306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
