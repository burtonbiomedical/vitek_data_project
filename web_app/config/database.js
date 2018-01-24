if(process.env.NODE_ENV === 'production'){
  module.exports = {mongoURL: 'mongodb://USER:PASSWORD@ds211588.mlab.com:11588/vitekdb-prod'}
}else{
  module.exports = {mongoURL: 'mongodb://localhost/vitekdb-dev'}
}